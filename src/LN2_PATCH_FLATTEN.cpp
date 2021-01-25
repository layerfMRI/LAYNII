#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PATCH_FLATTEN: Flatten a patch of cortex using 2D flat coordinate\n"
    "                   and cortical a depth measurement.\n"
    "\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "!!! WORK IN PROGRESS... EXPERIMENTAL PROGRAM. !!!\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "\n"
    "Usage:\n"
    "    LN2_PATCH_FLATTEN -input activation.nii -coord_uv uv_coord.nii -coord_d layers_equidist.nii -bins_u 50 -bins_v 50\n"
    "    LN2_PATCH_FLATTEN -input activation.nii -coord_uv uv_coord.nii -coord_d metric_equidist.nii -bins_u 50 -bins_v 50 -bins_d 21\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -input    : Nifti image with values that will be projected on a flat image.\n"
    "                For example an activation map or another measurement like curvature.\n"
    "    -coord_uv : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -coord_d  : A 3D nifti file that contains cortical depth measurements or layers.\n"
    "                For example either LN2_LAYERS output named 'layers' or 'metric'.\n"
    "    -domain   : A 3D binary nifti file to limit the flattened voxels.\n"
    "                For example LN2_MULTILATERATE output named 'perimeter_chunk.'\n"
    "    -bins_u   : Number of bins for the flat image U coordinate.\n"
    "    -bins_v   : Number of bins for the flat image V coordinate.\n"
    "    -bins_d   : (Optional) Number of bins for the flat image D coordinate.\n"
    "                Only use if '-coord_d' input is a metric file.\n"
    "    -debug    : (Optional) Save extra intermediate outputs.\n"
    "    -output   : (Optional) Output basename for all outputs.\n"
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
    int bins_u = 10, bins_v = 10, bins_d = 1;
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
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin4 = argv[ac];
        } else if (!strcmp(argv[ac], "-bins_u")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -bins_u\n");
                return 1;
            }
            bins_u = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-bins_v")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -bins_v\n");
                return 1;
            }
            bins_v = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-bins_d")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -bins_d\n");
                return 1;
            }
            bins_d = atof(argv[ac]);
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

    // ========================================================================
    // Determine the type of depth file
    float min_d = std::numeric_limits<float>::max();
    float max_d = std::numeric_limits<float>::min();

    // Check D coordinate min & max
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(coords_d_data + i) != 0) {
            if (*(coords_d_data + i) < min_d) {
                min_d = *(coords_d_data + i);
            }
            if (*(coords_d_data + i) > max_d) {
                max_d = *(coords_d_data + i);
            }
        }
    }

    // Determine whether depth input is a metric file or a layer file
    bool mode_depth_metric = NULL;
    if (min_d >= 0 && max_d <= 1) {
        cout << "  Depth input is a metric file (values are in between 0-1)." << endl;
        mode_depth_metric = true;
    } else if (min_d >= 0) {
        cout << "  Depth input is a layer file (values are positive integers)." << endl;
        mode_depth_metric = false;
    } else {
        cout << "  ERROR! Depth input contains negative values!" << endl;
        return 1;
    }

    // ========================================================================
    // Prepare outputs
    nifti_image* out_cells = copy_nifti_as_int32(nii_input);
    int32_t* out_cells_data = static_cast<int32_t*>(out_cells->data);

    for (int i = 0; i != nr_voxels; ++i) {
        *(out_cells_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // Determine flat image dimensions
    int nr_cells = bins_u * bins_v;
    if (mode_depth_metric == false) {  // Layer file
        bins_d = max_d;
    }
    int nr_bins = nr_cells * bins_d;

    // Allocating new nifti for flat images
    nifti_image* flat_cells = nifti_copy_nim_info(nii1);
    flat_cells->datatype = NIFTI_TYPE_INT32;
    flat_cells->dim[0] = 4;  // For proper 4D nifti
    flat_cells->dim[1] = bins_u;
    flat_cells->dim[2] = bins_v;
    flat_cells->dim[3] = bins_d;
    flat_cells->dim[4] = 1;
    nifti_update_dims_from_array(flat_cells);
    flat_cells->nvox = nr_bins;
    flat_cells->nbyper = sizeof(int32_t);
    flat_cells->data = calloc(flat_cells->nvox, flat_cells->nbyper);
    flat_cells->scl_slope = 1;
    int32_t* flat_cells_data = static_cast<int32_t*>(flat_cells->data);

    for (int i = 0; i != nr_bins; ++i) {
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
    cout << "  U coordinate min & max: " << min_u << " | " << max_u << endl;
    cout << "  V coordinate min & max: " << min_v << " | " << max_v << endl;

    // ========================================================================
    // Visit each voxel to check their coordinate
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);

        float u = *(coords_uv_data + nr_voxels*0 + i);
        float v = *(coords_uv_data + nr_voxels*1 + i);

        // Normalize coordinates to 0-1 range
        u = (u - min_u) / (max_u + std::numeric_limits<float>::min() - min_u);
        v = (v - min_v) / (max_v + std::numeric_limits<float>::min() - min_v);
        // Scale with grid size
        u *= static_cast<float>(bins_u);
        v *= static_cast<float>(bins_v);
        // Cast to integer (floor & cast)
        int cell_idx_u = static_cast<int>(u);
        int cell_idx_v = static_cast<int>(v);

        // Handle depth separately
        float d = static_cast<float>(*(coords_d_data + i));
        int cell_idx_d = 0;
        if (mode_depth_metric) {  // Metric file
            if (d >= 1) {  // Include 1 in the max index
                cell_idx_d = bins_d;
            } else {  // Scale up and floor
                d *= bins_d;
                cell_idx_d = static_cast<int>(d);
            }
        } else {  // Layer file
            cell_idx_d = static_cast<int>(d - 1);
        }

        // Flat image cell index
        int j = bins_u * cell_idx_v + cell_idx_u;
        int k = cell_idx_d * nr_cells + j;

        // Write cell index to output
        *(out_cells_data + i) = j;

        // Write visited voxel value to flat cell
        *(flat_input_data + k) += *(nii_input_data + i);
        *(flat_density_data + k) += 1;
    }

    save_output_nifti(fout, "UV_cells", out_cells, true);

    // Take the mean of each projected cell value
    for (int i = 0; i != nr_bins; ++i) {
        if (*(flat_density_data + i) > 1) {
            *(flat_input_data + i) /= *(flat_density_data + i);
        }
    }

    save_output_nifti(fout, "flat_values", flat_input, true);
    save_output_nifti(fout, "flat_density", flat_density, true);

    cout << "\n  Finished." << endl;
    return 0;
}
