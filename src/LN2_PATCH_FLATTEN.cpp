#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PATCH_FLATTEN: Flatten a patch of cortex using 2D flat coordinates (U and V)\n"
    "                   and cortical a depth measurement (D). Intended to be used in\n"
    "                   combination with LN2_MULTILATERATE, and LN2_LAYERS outputs.\n"
    "\n"
    "Usage:\n"
    "    LN2_PATCH_FLATTEN -values curvature.nii -coord_uv uv_coord.nii -coord_d metric_equidist.nii -domain perimeter_chunk.nii -bins_u 50 -bins_v 50 -bins_d 21\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -values    : Nifti image with values that will be projected onto flat image.\n"
    "                 For example an activation map or another measurement like curvature.\n"
    "                 If this file is 4D (e.g. time series), the output will also be 4D.\n"
    "    -coord_uv  : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                 For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -coord_d   : A 3D nifti file that contains cortical depth measurements or layers.\n"
    "                 For example either LN2_LAYERS output named 'layers' or 'metric'.\n"
    "    -domain    : A 3D binary nifti file to limit the flattened voxels.\n"
    "                 For example LN2_MULTILATERATE output named 'perimeter_chunk.'\n"
    "    -bins_u    : Number of bins for the flat image U coordinate.\n"
    "    -bins_v    : Number of bins for the flat image V coordinate.\n"
    "    -bins_d    : (Optional) Number of bins for the flat image D coordinate.\n"
    "                 Only use if '-coord_d' input is a metric file.\n"
    "    -voronoi   : (Optional) Fill empty bin in flat image using Voronoi propagation.\n"
    "                 Same as nearest neighbour filling in the empty bins.\n"
    "    -density   : (Optional) Additional output showing how many voxel fall into\n"
    "                 the same flat bin.\n"
    "    -norm_mask : (Optional) Mask out flat domain voxels using L2 norm of coordinates.\n"
    "    -debug     : (Optional) Save extra intermediate outputs.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n"
    "Citation:\n"
    "    - Gulban, O. F., Bollmann, S., Huber, R., Wagstyl, K., Goebel, R., Poser,\n"
    "      B. A., Kay, K., Ivanov, D. (2022). Mesoscopic in vivo human T2* dataset\n"
    "      acquired using quantitative MRI at 7 Tesla. Neuroimage.\n"
    "      <https://doi.org/10.1016/j.neuroimage.2022.119733>\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL, *nii4 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL, *fin4=NULL;
    int ac;
    int bins_u = 10, bins_v = 10, bins_d = 1;
    bool mode_debug = false, mode_voronoi = false, mode_norm_mask = false;
    bool mode_density = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
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
        } else if (!strcmp(argv[ac], "-voronoi")) {
            mode_voronoi = true;
        } else if (!strcmp(argv[ac], "-density")) {
            mode_density = true;
        } else if (!strcmp(argv[ac], "-norm_mask")) {
            mode_norm_mask = true;
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
        fprintf(stderr, "** missing option '-values'\n");
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
    const int64_t size_x = nii1->nx;
    const int64_t size_y = nii1->ny;
    const int64_t size_z = nii1->nz;
    const int64_t size_time = nii1->nt;
    const int64_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
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
    // ========================================================================
    float min_d = std::numeric_limits<float>::max();
    float max_d = std::numeric_limits<float>::min();

    // Check D coordinate min & max
    for (int64_t i = 0; i != nr_voxels; ++i) {
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
    bool mode_depth_metric = false;
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
    // ========================================================================
    nifti_image* out_cells = copy_nifti_as_int32(domain);
    int32_t* out_cells_data = static_cast<int32_t*>(out_cells->data);

    for (int64_t i = 0; i != nr_voxels; ++i) {
        *(out_cells_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // Determine flat image dimensions
    // ------------------------------------------------------------------------
    int nr_cells = bins_u * bins_v;
    if (mode_depth_metric == false) {  // Layer file
        bins_d = max_d;
    }
    int nr_bins = nr_cells * bins_d;

    // Add bin dimensions into the output tag
    std::ostringstream tag_u, tag_v, tag_d;
    tag_u << bins_u;
    tag_v << bins_v;
    tag_d << bins_d;

    // Allocating new 4D nifti for flat images
    nifti_image* flat_4D = nifti_copy_nim_info(nii1);
    flat_4D->datatype = NIFTI_TYPE_INT32;
    flat_4D->dim[0] = 4;  // For proper 4D nifti
    flat_4D->dim[1] = bins_u;
    flat_4D->dim[2] = bins_v;
    flat_4D->dim[3] = bins_d;
    flat_4D->dim[4] = size_time;
    flat_4D->pixdim[1] = 1;
    flat_4D->pixdim[2] = 1;
    flat_4D->pixdim[3] = 1;
    nifti_update_dims_from_array(flat_4D);
    flat_4D->nvox = nr_bins * size_time;
    flat_4D->nbyper = sizeof(int32_t);
    flat_4D->data = calloc(flat_4D->nvox, flat_4D->nbyper);
    flat_4D->scl_slope = 1;
    int32_t* flat_4D_data = static_cast<int32_t*>(flat_4D->data);

    for (int i = 0; i != nr_bins*size_time; ++i) {
        *(flat_4D_data + i) = 0;
    }

    // Flat input projection
    nifti_image* flat_values = copy_nifti_as_float32(flat_4D);
    float* flat_values_data = static_cast<float*>(flat_values->data);

    // ------------------------------------------------------------------------
    // Allocating new 3D nifti for flat images
    nifti_image* flat_3D = nifti_copy_nim_info(nii1);
    flat_3D->datatype = NIFTI_TYPE_INT32;
    flat_3D->dim[0] = 4;  // For proper 4D nifti
    flat_3D->dim[1] = bins_u;
    flat_3D->dim[2] = bins_v;
    flat_3D->dim[3] = bins_d;
    flat_3D->dim[4] = 1;
    flat_3D->pixdim[1] = 1;
    flat_3D->pixdim[2] = 1;
    flat_3D->pixdim[3] = 1;
    nifti_update_dims_from_array(flat_3D);
    flat_3D->nvox = nr_bins;
    flat_3D->nbyper = sizeof(int32_t);
    flat_3D->data = calloc(flat_3D->nvox, flat_3D->nbyper);
    flat_3D->scl_slope = 1;
    int32_t* flat_3D_data = static_cast<int32_t*>(flat_3D->data);

    for (int i = 0; i != nr_bins; ++i) {
        *(flat_3D_data + i) = 0;
    }

    // Flat projection density map
    nifti_image* flat_density = copy_nifti_as_float32(flat_3D);
    float* flat_density_data = static_cast<float*>(flat_density->data);
    // Flat domain
    nifti_image* flat_domain = copy_nifti_as_float32(flat_3D);
    float* flat_domain_data = static_cast<float*>(flat_domain->data);

    // ------------------------------------------------------------------------
    // Allocating new 4D nifti for saveing the folded image coordinates in the 
    // flat image format. This is for back projection from flat to folded.
    nifti_image* flat_coords = nifti_copy_nim_info(nii1);
    flat_coords->datatype = NIFTI_TYPE_FLOAT32;
    flat_coords->dim[0] = 4;  // For proper 4D nifti
    flat_coords->dim[1] = bins_u;
    flat_coords->dim[2] = bins_v;
    flat_coords->dim[3] = bins_d;
    flat_coords->dim[4] = 3;
    flat_coords->pixdim[1] = 1;
    flat_coords->pixdim[2] = 1;
    flat_coords->pixdim[3] = 1;
    nifti_update_dims_from_array(flat_coords);
    flat_coords->nvox = nr_bins * 3;
    flat_coords->nbyper = sizeof(float);
    flat_coords->data = calloc(flat_coords->nvox, flat_coords->nbyper);
    flat_coords->scl_slope = 1;
    float* flat_coords_data = static_cast<float*>(flat_coords->data);

    for (int i = 0; i != nr_bins*3; ++i) {
        *(flat_coords_data + i) = 0;
    }

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
    // ========================================================================
    float min_u = std::numeric_limits<float>::max();
    float max_u = std::numeric_limits<float>::min();
    float min_v = std::numeric_limits<float>::max();
    float max_v = std::numeric_limits<float>::min();

    for (int64_t ii = 0; ii != nr_voi; ++ii) {
        int64_t i = *(voi_id + ii);
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
    // ========================================================================
    for (int64_t t = 0; t != size_time; ++t) {
        for (int64_t ii = 0; ii != nr_voi; ++ii) {
            int64_t i = *(voi_id + ii);

            float u = *(coords_uv_data + nr_voxels*0 + i);
            float v = *(coords_uv_data + nr_voxels*1 + i);

            // Normalize coordinates to 0-1 range
            u = (u - min_u) / (max_u + std::numeric_limits<float>::min() - min_u);
            v = (v - min_v) / (max_v + std::numeric_limits<float>::min() - min_v);
            // Scale with grid size
            u *= static_cast<float>(bins_u);
            v *= static_cast<float>(bins_v);
            // Cast to integer (floor & cast)
            int64_t cell_idx_u = static_cast<int64_t>(u);
            int64_t cell_idx_v = static_cast<int64_t>(v);

            // Handle depth separately
            float d = static_cast<float>(*(coords_d_data + i));
            int64_t cell_idx_d = 0;
            if (mode_depth_metric) {  // Metric file
                if (d >= 1) {  // Include 1 in the max index
                    cell_idx_d = bins_d;
                } else {  // Scale up and floor
                    d *= bins_d;
                    cell_idx_d = static_cast<int64_t>(d);
                }
            } else {  // Layer file
                cell_idx_d = static_cast<int64_t>(d - 1);
            }

            // Flat image cell index
            int64_t j = bins_u * cell_idx_v + cell_idx_u;
            int64_t k = cell_idx_d * nr_cells + j;

            // Write cell index to output
            *(out_cells_data + i) = j + 1;

            // Write visited voxel value to flat cell
            *(flat_values_data + k + t*nr_bins) += *(nii_input_data + i + t*nr_voxels);

            // Project folded data coordinates
            int64_t ix, iy, iz;
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
            *(flat_coords_data + k + nr_bins*0) += static_cast<float>(ix);
            *(flat_coords_data + k + nr_bins*1) += static_cast<float>(iy);
            *(flat_coords_data + k + nr_bins*2) += static_cast<float>(iz);

            if (t==0) {  // Write 3D values once
                *(flat_density_data + k) += 1;
                *(flat_domain_data + k) += *(domain_data + i);
            }
        }
    }

    // Take the mean of each projected cell value
    for (int64_t t = 0; t != size_time; ++t) {
        for (int64_t i = 0; i != nr_bins; ++i) {
            if (*(flat_density_data + i) > 1) {
                *(flat_values_data + i + t*nr_bins) /= *(flat_density_data + i);

                *(flat_coords_data + i + nr_bins*0) /= *(flat_density_data + i);
                *(flat_coords_data + i + nr_bins*1) /= *(flat_density_data + i);
                *(flat_coords_data + i + nr_bins*2) /= *(flat_density_data + i);

                if (t==0) {  // Do once for 3D values
                    *(flat_domain_data + i) /= *(flat_density_data + i);
                    // Ceil domain average to ensure the edges are prioritized
                    *(flat_domain_data + i) = std::ceil(*(flat_domain_data + i));
                }
            }
        }
    }

    // ========================================================================
    // Optional Voronoi filling for empty flat bins
    // ========================================================================
    if (mode_voronoi) {
        cout << "\n  Start Voronoi (nearest neighbor) filling-in..." << endl;

        // Prepare additional flat niftis
        nifti_image* flood_step = copy_nifti_as_float32(flat_4D);
        float* flood_step_data = static_cast<float*>(flood_step->data);
        nifti_image* flood_dist = copy_nifti_as_float32(flat_4D);
        float* flood_dist_data = static_cast<float*>(flood_dist->data);

        // --------------------------------------------------------------------
        const int64_t size_x = bins_u;
        const int64_t size_y = bins_v;
        const int64_t size_z = bins_d;
        const int64_t end_x = size_x - 1;
        const int64_t end_y = size_y - 1;
        const int64_t end_z = size_z - 1;

        const float dX = 1;
        const float dY = 1;
        const float dZ = 1;

        // Short diagonals
        const float dia_xy = sqrt(dX * dX + dY * dY);
        const float dia_xz = sqrt(dX * dX + dZ * dZ);
        const float dia_yz = sqrt(dY * dY + dZ * dZ);
        // Long diagonals
        const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

        for (int64_t t=0; t!=size_time; ++t) {
            // Initialize grow volume
            for (int64_t i = 0; i != nr_bins; ++i) {
                if (*(flat_values_data + i + t*nr_bins) != 0) {
                    *(flood_step_data + i) = 1.;
                    *(flood_dist_data + i) = 0.;
                } else {
                    *(flood_step_data + i) = 0.;
                    *(flood_dist_data + i) = 0.;
                }
            }
            int64_t grow_step = 1, bin_counter = 1;
            int64_t ix, iy, iz, j;
            float d;
            if (size_time > 1) {
                cout << "  Doing 4th dimension: " << t+1 << endl;
            }
            bin_counter = nr_bins;

            while (bin_counter != 0) {
                bin_counter = 0;
                for (int64_t i = 0; i != nr_bins; ++i) {
                    if (*(flood_step_data + i) == grow_step) {
                        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                        bin_counter += 1;

                        // --------------------------------------------------------
                        // 1-jump neighbours
                        // --------------------------------------------------------
                        if (ix > 0) {
                            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x) {
                            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy > 0) {
                            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy < end_y) {
                            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iz > 0) {
                            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iz < end_z) {
                            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }

                        // --------------------------------------------------------
                        // 2-jump neighbours
                        // --------------------------------------------------------

                        if (ix > 0 && iy > 0) {
                            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iy < end_y) {
                            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy > 0) {
                            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy < end_y) {
                            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iz > 0) {
                            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iz > 0) {
                            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }

                        // --------------------------------------------------------
                        // 3-jump neighbours
                        // --------------------------------------------------------
                        if (ix > 0 && iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix > 0 && iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                        if (ix < end_x && iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(flat_density_data + j) = *(flat_density_data + i);
                                *(flat_domain_data + j) = *(flat_domain_data + i);
                                *(flat_coords_data + nr_bins*0 + j) = *(flat_coords_data + nr_bins*0 + i);
                                *(flat_coords_data + nr_bins*1 + j) = *(flat_coords_data + nr_bins*1 + i);
                                *(flat_coords_data + nr_bins*2 + j) = *(flat_coords_data + nr_bins*2 + i);
                            }
                        }
                    }
                }
            grow_step += 1;
            }
        }

        // NOTE(Option 2) Mask values based on radius
        if (mode_norm_mask) {
            for (int64_t i = 0; i != nr_bins; ++i) {
                float coord_u = i % bins_u;
                float coord_v = floor(i % nr_cells / bins_v);
                coord_u /= bins_u;
                coord_v /= bins_v;
                coord_u -= 0.5;
                coord_v -= 0.5;

                float mag_uv = sqrt(pow(coord_u, 2) + pow(coord_v, 2));
                if (mag_uv > 0.5) {
                    *(flat_domain_data + i) = 2;
                }
            }
        }

        // NOTE(Option 1): Mask values outside of the flattened disk
        for (int64_t i = 0; i != nr_bins; ++i) {
            if ( *(flat_domain_data + i) != 1 ) {
                *(flat_density_data + i) = 0;
                *(flat_coords_data + nr_bins*0 + i) = 0;
                *(flat_coords_data + nr_bins*1 + i) = 0;
                *(flat_coords_data + nr_bins*2 + i) = 0;
            }
        }

        for (int64_t i = 0; i != nr_bins; ++i) {
            for (int64_t t = 0; t != size_time; ++t) {
                if ( *(flat_domain_data + i) != 1 ) {
                    *(flat_values_data + nr_bins*t + i) = 0;
                }
            }
        }

        save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_voronoi", flat_values, true);
        save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_foldedcoords_voronoi", flat_coords, true);
        if (mode_density) {
            save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_density_voronoi", flat_density, true);
        }
        if (mode_debug) {
            save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_domain_voronoi", flat_domain, true);
        }
    } else {
        if (mode_debug) {
            save_output_nifti(fout, "UV_bins_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str(), out_cells, true);
        }
        save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str(), flat_values, true);
        save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_foldedcoords", flat_coords, true);
        if (mode_density) {
            save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_density", flat_density, true);
        }
        if (mode_debug) {
            save_output_nifti(fout, "flat_"+tag_u.str()+"x"+tag_v.str()+"x"+tag_d.str()+"_domain", flat_domain, true);
        }

    }

    cout << "\n  Finished." << endl;
    return 0;
}
