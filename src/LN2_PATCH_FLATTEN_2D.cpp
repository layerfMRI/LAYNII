#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PATCH_FLATTEN_2D: Equivalent of LN2_PATCH_FLATTEN but designed to work with 2D\n"
    "                      images such as histology images. Flattens a patch of cortex two\n"
    "                      coordinates. For instance, one coordinate indicates the radial direction\n"
    "                      with regards to the cortical surface (e.g. 'cortical depths' or 'layers'.\n"
    "                      the other coordinate indicates the tangential direction with regards to\n"
    "                      the cortical surface. Radial coordinate can be the 'metric' or 'layers'\n"
    "                      output of LN2_LAYERS. Tangential coordinate can be the output of \n"
    "                      e.g. LN2_GEODISTANCE.\n"
    "\n"
    "Usage:\n"
    "    LN2_PATCH_FLATTEN_2D -values curvature.nii -coord_tan dist_columnar.nii -coord_rad metric_equidist.nii -domain mask.nii -bins_rad 21 -bins_tan 50\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -values    : Nifti image with values that will be projected onto flat image.\n"
    "                 For example an activation map or another measurement like curvature.\n"
    "                 If this file is 4D (e.g. time series), the output will also be 4D.\n"
    "    -coord_tan : A 2D nifti file that contains *tangential* (w.r.t. cortical surface)\n"
    "                 distance measurements. For example, LN2_GEODISTANCE output distances.\n"
    "    -coord_rad : A 2D nifti file that contains *radial* (w.r.t. cortical surface).\n"
    "                 distance measurements. For example, LN2_LAYERS output named 'metric'.\n"
    "    -domain    : A 2D binary nifti file to limit the flattened voxels.\n"
    "    -bins_rad  : Number of bins for the flat image radial coordinate.\n"
    "    -bins_tan  : Number of bins for the flat image tangential coordinate.\n"
    "    -voronoi   : (Optional) Fill empty bin in flat image using Voronoi propagation.\n"
    "                 Same as nearest neighbour filling in the empty bins.\n"
    "    -density   : (Optional) Additional output showing how many voxel fall into\n"
    "                 the same flat bin.\n"
    "    -debug     : (Optional) Save extra intermediate outputs.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL, *nii4 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL, *fin4=NULL;
    int ac;
    int bins_rad = 21, bins_tan = 50;
    bool mode_debug = false, mode_voronoi = false, mode_density = false;

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
        } else if (!strcmp(argv[ac], "-coord_tan")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_tan\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-coord_rad")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_rad\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin4 = argv[ac];
        } else if (!strcmp(argv[ac], "-bins_rad")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -bins_rad\n");
                return 1;
            }
            bins_rad = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-bins_tan")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -bins_tan\n");
                return 1;
            }
            bins_tan = atof(argv[ac]);
         } else if (!strcmp(argv[ac], "-voronoi")) {
            mode_voronoi = true;
        } else if (!strcmp(argv[ac], "-density")) {
            mode_density = true;
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
        fprintf(stderr, "** missing option '-coords_tangential'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-coords_radial'\n");
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
    const int size_time = nii1->nt;
    const int nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* coords_tan = copy_nifti_as_float32(nii2);
    float* coords_tan_data = static_cast<float*>(coords_tan->data);
    nifti_image* coords_rad = copy_nifti_as_float32(nii3);
    float* coords_rad_data = static_cast<float*>(coords_rad->data);
    nifti_image* domain = copy_nifti_as_int32(nii4);
    int32_t* domain_data = static_cast<int32_t*>(domain->data);

    // ========================================================================
    // Determine the type of depth file
    // ========================================================================
    float min_rad = std::numeric_limits<float>::max();
    float max_rad = std::numeric_limits<float>::min();

    // Check D coordinate min & max
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(coords_rad_data + i) != 0) {
            if (*(coords_rad_data + i) < min_rad) {
                min_rad = *(coords_rad_data + i);
            }
            if (*(coords_rad_data + i) > max_rad) {
                max_rad = *(coords_rad_data + i);
            }
        }
    }

    // Determine whether depth input is a metric file or a layer file
    bool mode_depth_metric = false;
    if (min_rad >= 0 && max_rad <= 1) {
        cout << "  Depth input is a metric file (values are in between 0-1)." << endl;
        mode_depth_metric = true;
    } else if (min_rad >= 0) {
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

    for (int i = 0; i != nr_voxels; ++i) {
        *(out_cells_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // Determine flat image dimensions
    // ------------------------------------------------------------------------
    if (mode_depth_metric == false) {  // Layer file
        bins_rad = max_rad;
    }
    int nr_bins = bins_tan * bins_rad;

    // Add bin dimensions into the output tag
    std::ostringstream tag_rad, tag_tan;
    tag_rad << bins_rad;
    tag_tan << bins_tan;

    // Allocating new 4D nifti for flat images
    nifti_image* flat_3D = nifti_copy_nim_info(nii1);
    flat_3D->datatype = NIFTI_TYPE_INT32;
    flat_3D->dim[0] = 4;  // For proper 4D nifti
    flat_3D->dim[1] = bins_tan;
    flat_3D->dim[2] = bins_rad;
    flat_3D->dim[3] = 1;
    flat_3D->dim[4] = size_time;
    flat_3D->pixdim[1] = 1;
    flat_3D->pixdim[2] = 1;
    flat_3D->pixdim[3] = 1;
    nifti_update_dims_from_array(flat_3D);
    flat_3D->nvox = nr_bins * size_time;
    flat_3D->nbyper = sizeof(int32_t);
    flat_3D->data = calloc(flat_3D->nvox, flat_3D->nbyper);
    flat_3D->scl_slope = 1;
    int32_t* flat_3D_data = static_cast<int32_t*>(flat_3D->data);

    for (int i = 0; i != nr_bins*size_time; ++i) {
        *(flat_3D_data + i) = 0;
    }

    // Flat input projection
    nifti_image* flat_values = copy_nifti_as_float32(flat_3D);
    float* flat_values_data = static_cast<float*>(flat_values->data);

    // ------------------------------------------------------------------------
    // Allocating new 3D nifti for flat images
    nifti_image* flat_2D = nifti_copy_nim_info(nii1);
    flat_2D->datatype = NIFTI_TYPE_INT32;
    flat_2D->dim[0] = 4;  // For proper 4D nifti
    flat_2D->dim[1] = bins_tan;
    flat_2D->dim[2] = bins_rad;
    flat_2D->dim[3] = 1;
    flat_2D->dim[4] = 1;
    flat_2D->pixdim[1] = 1;
    flat_2D->pixdim[2] = 1;
    flat_2D->pixdim[3] = 1;
    nifti_update_dims_from_array(flat_2D);
    flat_2D->nvox = nr_bins;
    flat_2D->nbyper = sizeof(int32_t);
    flat_2D->data = calloc(flat_2D->nvox, flat_2D->nbyper);
    flat_2D->scl_slope = 1;
    int32_t* flat_2D_data = static_cast<int32_t*>(flat_2D->data);

    for (int i = 0; i != nr_bins; ++i) {
        *(flat_2D_data + i) = 0;
    }

    // Flat projection density map
    nifti_image* flat_density = copy_nifti_as_float32(flat_2D);
    float* flat_density_data = static_cast<float*>(flat_density->data);
    // Flat domain
    nifti_image* flat_domain = copy_nifti_as_float32(flat_2D);
    float* flat_domain_data = static_cast<float*>(flat_domain->data);

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
    float min_tan = std::numeric_limits<float>::max();
    float max_tan = std::numeric_limits<float>::min();

    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        // Check tangential coordinate min & max
        if (*(coords_tan_data + i) < min_tan) {
            min_tan = *(coords_tan_data + i);
        }
        if (*(coords_tan_data + i) > max_tan) {
            max_tan = *(coords_tan_data + i);
        }
    }
    cout << "  Tangential coordinate min & max: " << min_tan << " | " << max_tan << endl;

    // ========================================================================
    // Visit each voxel to check their coordinate
    // ========================================================================
    for (int t = 0; t != size_time; ++t) {
        for (int ii = 0; ii != nr_voi; ++ii) {
            int i = *(voi_id + ii);

            float u = *(coords_tan_data + i);
            cout << "U: " << u << endl;

            // Normalize coordinates to 0-1 range
            u = (u - min_tan) / (max_tan + std::numeric_limits<float>::min() - min_tan);
            // Scale with grid size
            u *= static_cast<float>(bins_tan);
            // Cast to integer (floor & cast)
            int cell_idx_u = static_cast<int>(u);

            cout << "CELL_IDX_U: " << cell_idx_u << endl;

            // Handle depth separately
            float d = static_cast<float>(*(coords_rad_data + i));
            // cout << "D: " << d << endl;
            int cell_idx_d = 0;
            if (mode_depth_metric) {  // Metric file
                if (d >= 1) {  // Include 1 in the max index
                    cell_idx_d = bins_rad;
                } else {  // Scale up and floor
                    d *= bins_rad;
                    cell_idx_d = static_cast<int>(d);
                }
            } else {  // Layer file
                cell_idx_d = static_cast<int>(d - 1);
            }

            // cout << "CELL_IDX_D: " << cell_idx_d << endl;

            // Flat image cell index
            int j = bins_tan * cell_idx_d + cell_idx_u;

            // Write cell index to output
            *(out_cells_data + i) = j + 1;

            // Write visited voxel value to flat cell
            *(flat_values_data + j + t*nr_bins) += *(nii_input_data + i + t*nr_voxels);

            if (t==0) {  // Write 3D values once
                *(flat_density_data + j) += 1;
                *(flat_domain_data + j) += *(domain_data + i);
            }
        }
    }

    // Take the mean of each projected cell value
    for (int t = 0; t != size_time; ++t) {
        for (int i = 0; i != nr_bins; ++i) {
            if (*(flat_density_data + i) > 1) {
                *(flat_values_data + i + t*nr_bins) /= *(flat_density_data + i);

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

        // // Prepare additional flat niftis
        // nifti_image* flood_step = copy_nifti_as_float32(flat_3D);
        // float* flood_step_data = static_cast<float*>(flood_step->data);
        // nifti_image* flood_dist = copy_nifti_as_float32(flat_3D);
        // float* flood_dist_data = static_cast<float*>(flood_dist->data);

        // // --------------------------------------------------------------------
        // const int size_x = bins_rad;
        // const int size_y = bins_tan;
        // const int size_z = bins_rad;
        // const int end_x = size_x - 1;
        // const int end_y = size_y - 1;
        // const int end_z = size_z - 1;

        // const float dX = 1;
        // const float dY = 1;
        // const float dZ = 1;

        // // Short diagonals
        // const float dia_xy = sqrt(dX * dX + dY * dY);
        // const float dia_xz = sqrt(dX * dX + dZ * dZ);
        // const float dia_yz = sqrt(dY * dY + dZ * dZ);
        // // Long diagonals
        // const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

        // for (int t=0; t!=size_time; ++t) {
        //     // Initialize grow volume
        //     for (int i = 0; i != nr_bins; ++i) {
        //         if (*(flat_values_data + i + t*nr_bins) != 0) {
        //             *(flood_step_data + i) = 1.;
        //             *(flood_dist_data + i) = 0.;
        //         } else {
        //             *(flood_step_data + i) = 0.;
        //             *(flood_dist_data + i) = 0.;
        //         }
        //     }
        //     int grow_step = 1, bin_counter = 1;
        //     int ix, iy, iz, j;
        //     float d;
        //     if (size_time > 1) {
        //         cout << "  Doing 4th dimension: " << t+1 << endl;
        //     }
        //     bin_counter = nr_bins;

        //     while (bin_counter != 0) {
        //         bin_counter = 0;
        //         for (int i = 0; i != nr_bins; ++i) {
        //             if (*(flood_step_data + i) == grow_step) {
        //                 tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        //                 bin_counter += 1;

        //                 // --------------------------------------------------------
        //                 // 1-jump neighbours
        //                 // --------------------------------------------------------
        //                 if (ix > 0) {
        //                     j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dX;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x) {
        //                     j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dX;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy > 0) {
        //                     j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dY;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy < end_y) {
        //                     j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dY;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iz > 0) {
        //                     j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dZ;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iz < end_z) {
        //                     j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dZ;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }

        //                 // --------------------------------------------------------
        //                 // 2-jump neighbours
        //                 // --------------------------------------------------------

        //                 if (ix > 0 && iy > 0) {
        //                     j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xy;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iy < end_y) {
        //                     j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xy;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy > 0) {
        //                     j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xy;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy < end_y) {
        //                     j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xy;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy > 0 && iz > 0) {
        //                     j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_yz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy > 0 && iz < end_z) {
        //                     j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_yz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy < end_y && iz > 0) {
        //                     j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_yz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (iy < end_y && iz < end_z) {
        //                     j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_yz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iz > 0) {
        //                     j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iz > 0) {
        //                     j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iz < end_z) {
        //                     j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iz < end_z) {
        //                     j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }

        //                 // --------------------------------------------------------
        //                 // 3-jump neighbours
        //                 // --------------------------------------------------------
        //                 if (ix > 0 && iy > 0 && iz > 0) {
        //                     j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iy > 0 && iz < end_z) {
        //                     j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iy < end_y && iz > 0) {
        //                     j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy > 0 && iz > 0) {
        //                     j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix > 0 && iy < end_y && iz < end_z) {
        //                     j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy > 0 && iz < end_z) {
        //                     j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy < end_y && iz > 0) {
        //                     j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //                 if (ix < end_x && iy < end_y && iz < end_z) {
        //                     j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
        //                     d = *(flood_dist_data + i) + dia_xyz;
        //                     if (d < *(flood_dist_data + j)
        //                         || *(flood_dist_data + j) == 0) {
        //                         *(flat_values_data + j + t*nr_bins) = *(flat_values_data + i + t*nr_bins);
        //                         *(flood_dist_data + j) = d;
        //                         *(flood_step_data + j) = grow_step + 1;
        //                         *(flat_density_data + j) = *(flat_density_data + i);
        //                         *(flat_domain_data + j) = *(flat_domain_data + i);
        //                     }
        //                 }
        //             }
        //         }
        //     grow_step += 1;
        //     }
        // }

        // // NOTE(Option 2) Mask values based on radius
        // if (mode_norm_mask) {
        //     for (int i = 0; i != nr_bins; ++i) {
        //         float coord_u = i % bins_rad;
        //         float coord_v = floor(i % nr_cells / bins_tan);
        //         coord_u /= bins_rad;
        //         coord_v /= bins_tan;
        //         coord_u -= 0.5;
        //         coord_v -= 0.5;

        //         float mag_uv = sqrt(pow(coord_u, 2) + pow(coord_v, 2));
        //         if (mag_uv > 0.5) {
        //             *(flat_domain_data + i) = 2;
        //         }
        //     }
        // }

        // // NOTE(Option 1): Mask values outside of the flattened disk
        // for (int t = 0; t != size_time; ++t) {
        //     for (int i = 0; i != nr_bins; ++i) {
        //         if (*(flat_domain_data + i) != 1) {
        //             *(flat_values_data + i + t*nr_bins) = 0;
        //             *(flat_density_data + i) = 0;
        //             *(flat_coords_radata + nr_bins*0 + i) = 0;
        //             *(flat_coords_radata + nr_bins*1 + i) = 0;
        //             *(flat_coords_radata + nr_bins*2 + i) = 0;
        //         }
        //     }
        // }

        // save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"x"+tag_d.str()+"_voronoi", flat_values, true);
        // save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"x"+tag_d.str()+"_foldedcoords_voronoi", flat_coords, true);
        // if (mode_density) {
        //     save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"x"+tag_d.str()+"_density_voronoi", flat_density, true);
        // }
        // if (mode_debug) {
        //     save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"x"+tag_d.str()+"_domain_voronoi", flat_domain, true);
        // }
    } else {
        save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str(), flat_values, true);
        if (mode_density) {
            save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"_density", flat_density, true);
        }
        if (mode_debug) {
            save_output_nifti(fout, "flat_"+tag_rad.str()+"x"+tag_tan.str()+"_domain", flat_domain, true);
        }

    }

    cout << "\n  Finished." << endl;
    return 0;
}
