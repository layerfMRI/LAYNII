#include "../dep/laynii_lib.h"
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>

int show_help(void) {
    printf(
    "LN2_UVD_FILTER: Filter using flat coordinates (UV) and depth (D). This program:\n"
    "                  1. Evaluates voxels within a cylindrical window centered\n"
    "                     at each UVD coordinate.\n"
    "                  2. Performs the chosen operation (median, min, max, ...)\n"
    "                     within each window.\n"
    "                  3. Writes the result back to the center voxel of each window.\n"
    "\n"
    "Usage:\n"
    "    LN2_UVD_FILTER -values activation.nii -coord_uv uv_coord.nii -coord_d metric_equidist.nii -domain mask.nii -radius 3 -height 0.25\n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -values        : Nifti image with values that will be filtered.\n"
    "                     For example an activation map or anatomical T1w images.\n"
    "    -coord_uv      : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                     For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -coord_d       : A 3D nifti file that contains cortical depth measurements or layers.\n"
    "                     For example either LN2_LAYERS output named 'metric'.\n"
    "    -domain        : A 3D binary nifti file to limit the flattened voxels.\n"
    "                     For example LN2_MULTILATERATE output named 'perimeter_chunk.'\n"
    "    -radius        : Radius of cylinder that will be passed over UV coordinates.\n"
    "                     In units of UV coordinates, which often are in mm.\n"
    "    -height        : height/height of cylinder that will be passed over D (depth)\n"
    "                     coordinates. In units of normalized depth metric, which\n"
    "                     are often in 0-1 range. If you intend to make a cylinder\n"
    "                     that covers all depths at all times, make sure to enter '2'\n"
    "                     as the cylinder height.\n"
    "    -median        : (Default) Take the median within the window.\n"
    "    -min           : (Optional) Take the minimum within the window.\n"
    "    -max           : (Optional) Take the maximum within the window.\n"
    "    -columns       : (Optional) Take the mode within the window.\n"
    "    -peak_d        : (Optional) Take depth of the maximum value in the window.\n"
    "    -count_uniques : (Optional) Count number of uniquely labeled voxels within\n"
    "                     each window.\n"
    "    -output        : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL, *nii4 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL, *fin4=NULL;
    int ac;
    float radius = 3, height = 0.25;
    bool mode_median = true, mode_min = false, mode_max = false;
    bool mode_cols = false, mode_peak = false, mode_count_uniques = false;

    // Process user options
    if (argc < 2) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strcmp(argv[ac], "-help")) {
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
        } else if (!strcmp(argv[ac], "-radius")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radius\n");
                return 1;
            }
            radius = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-height")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -height\n");
                return 1;
            }
            height = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-median")) {
            mode_median = true;
            mode_max = false;
            mode_cols = false;
            mode_min = false;
            mode_peak = false;
            mode_count_uniques = false;
        } else if (!strcmp(argv[ac], "-min")) {
            mode_median = false;
            mode_max = false;
            mode_cols = false;
            mode_min = true;
            mode_peak = false;
            mode_count_uniques = false;
        } else if (!strcmp(argv[ac], "-max")) {
            mode_median = false;
            mode_min = false;
            mode_cols = false;
            mode_max = true;
            mode_peak = false;
            mode_count_uniques = false;
        } else if (!strcmp(argv[ac], "-columns")) {
            mode_median = false;
            mode_min = false;
            mode_max = false;
            mode_cols = true;
            mode_peak = false;
            mode_count_uniques = false;
        } else if (!strcmp(argv[ac], "-peak_d")) {
            mode_median = false;
            mode_min = false;
            mode_cols = false;
            mode_max = false;
            mode_peak = true;
            mode_count_uniques = false;
        } else if (!strcmp(argv[ac], "-count_uniques")) {
            mode_median = false;
            mode_min = false;
            mode_cols = false;
            mode_max = false;
            mode_peak = false;
            mode_count_uniques = true;
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
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin4);
        return 2;
    }

    log_welcome("LN2_UVD_FILTER");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);
    log_nifti_descriptives(nii4);

    // Get dimensions of input
    const int nr_voxels = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* coords_uv = copy_nifti_as_float32(nii2);
    float* coords_uv_data = static_cast<float*>(coords_uv->data);
    nifti_image* coords_d = copy_nifti_as_float32(nii3);
    float* coords_d_data = static_cast<float*>(coords_d->data);
    nifti_image* domain = copy_nifti_as_float32(nii4);
    float* domain_data = static_cast<float*>(domain->data);

    // ========================================================================
    // Prepare outputs
    nifti_image* nii_output = copy_nifti_as_float32(nii_input);
    float* nii_output_data = static_cast<float*>(nii_output->data);

    // For extra output
    nifti_image* nii_output_extra = copy_nifti_as_float32(nii_input);
    float* nii_output_extra_data = static_cast<float*>(nii_output_extra->data);

    nifti_image* temp_nii_output_extra = copy_nifti_as_float32(nii_input);
    float* temp_nii_output_extra_data = static_cast<float*>(temp_nii_output_extra->data);

    if (mode_min) {
        for (int i = 0; i != nr_voxels; ++i) {
            *(nii_output_data + i) = 0;
        }
    }


    if (mode_cols || mode_max) {
        for (int i = 0; i != nr_voxels; ++i) {
            *(nii_output_data + i) = 0;
            *(nii_output_extra_data + i) = 0;
            *(temp_nii_output_extra_data + i) = 0;
        }
    }

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    int nr_voi = 0;  // Voxels of interest
    vector <int> vec_voi_id;
    vector <float> vec_u, vec_v, vec_d, vec_val;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(domain_data + i) != 0){
            vec_voi_id.push_back(i);
            vec_u.push_back(*(coords_uv_data + nr_voxels*0 + i));
            vec_v.push_back(*(coords_uv_data + nr_voxels*1 + i));
            vec_d.push_back(*(coords_d_data + i));
            vec_val.push_back(*(nii_input_data + i));
            nr_voi += 1;
        }
    }

    // ========================================================================
    // Visit each voxel to check their coordinate
    // ========================================================================
    float half_height = height / 2;
    float radius_sqr = radius * radius;
    for (int i = 0; i != nr_voi; ++i) {
        cout << "\r    " << i * 100 / nr_voi << " %" << flush;
        vector <float> temp_vec;
        vector <int> temp_vec_id;
        vector <float> temp_vec_d;

        // --------------------------------------------------------------------
        // Cylinder windowing in UVD space
        // --------------------------------------------------------------------
        for (int j = 0; j != nr_voi; ++j) {
            // Compute distances relative to reference UVD
            if (abs(vec_d[i] - vec_d[j]) < half_height) {  // Check height
                float dist_uv = (vec_u[i] - vec_u[j])*(vec_u[i] - vec_u[j])
                    + (vec_v[i] - vec_v[j])*(vec_v[i] - vec_v[j]);
                if (dist_uv < radius_sqr) {  // Check Euclidean distance
                    temp_vec.push_back(vec_val[j]);
                    temp_vec_id.push_back(vec_voi_id[j]);
                    temp_vec_d.push_back(vec_d[j]);

                }
            }
        }

        int n = temp_vec.size();
        // --------------------------------------------------------------------
        // Find median
        // --------------------------------------------------------------------
        if (mode_median) {
            float m;
            if (n % 2 == 0) {  // even
                std::nth_element(temp_vec.begin(),
                temp_vec.begin() + n / 2,
                temp_vec.end());

                std::nth_element(temp_vec.begin(),
                temp_vec.begin() + (n - 1) / 2,
                temp_vec.end());

                m = (temp_vec[n / 2] + temp_vec[(n - 1) / 2]) / 2.0;

            } else {  // odd
                std::nth_element(temp_vec.begin(),
                temp_vec.begin() + n / 2,
                temp_vec.end());

                m = temp_vec[n / 2];
            }

            *(nii_output_data + vec_voi_id[i]) = m;
        }

        // --------------------------------------------------------------------
        // Find minimum
        // --------------------------------------------------------------------
        if (mode_min) {
            float temp_ref = vec_val[i];
            float temp_min = vec_val[i];
            for (int j = 0; j != n; ++j) {
                if (temp_vec[j] < temp_min) {
                    temp_min = temp_vec[j];
                }
            }

            if (temp_min < temp_ref) {
                *(nii_output_data + vec_voi_id[i]) = 0;
            } else {
                *(nii_output_data + vec_voi_id[i]) = 1;
            }
        }

        // --------------------------------------------------------------------
        // Find maximum (Works with binary mask)
        // NOTE: temp_mask = 1 or 0
        // --------------------------------------------------------------------
        if (mode_max) {
            float temp_ref = vec_val[i];
            float temp_max = vec_val[i];
            for (int j = 0; j != n; ++j) {
                if (temp_vec[j] > temp_max) {
                    temp_max = temp_vec[j];
                }
            }
            *(nii_output_data + vec_voi_id[i]) = temp_max;
            *(temp_nii_output_extra_data + vec_voi_id[i]) = static_cast<float>(n);
        }

        // --------------------------------------------------------------------
        // A) Find functional columns: write back to a single voxel
        // TODO[Faruk]: Code review this part.
        // --------------------------------------------------------------------
        if (mode_cols) {
            int count1 = 0, count2 = 0,  count3 = 0, count4 = 0;
            int m, c, t;
            // Count occurrences
            for (int j = 0; j != n; ++j) {
                if (temp_vec[j] == 1) {
                    count1 += 1;
                } else if (temp_vec[j] == 2) {
                    count2 += 1;
                } else if (temp_vec[j] == 3) {
                    count3 += 1;
                } else if (temp_vec[j] == 4) {
                    count4 += 1;
                }
            }

            // Find the maximum (most common label)
            if (count1 > count2 && count1 > count3 && count1 > count4) {
                m = 1; c = count1;
            } else if (count2 > count1 && count2 > count3 && count2 > count4) {
                m = 2; c = count2;
            } else if (count3 > count1 && count3 > count2 && count3 > count4) {
                m = 3; c = count3;
            } else if (count4 > count1 && count4 > count2 && count4 > count3) {
                m = 4; c = count4;
            }
            t = count1 + count2 + count3 + count4;

            // Write the output (voxel wise)
            *(nii_output_data + vec_voi_id[i]) = m;
            *(nii_output_extra_data + vec_voi_id[i]) = static_cast<float>(c)/t;
            *(temp_nii_output_extra_data + vec_voi_id[i]) = static_cast<float>(t);
        }
        // --------------------------------------------------------------------
        // B) Find functional columns (strict definition): Write back to all
        // window
        // TODO[Faruk]: This section might be redundant or needs to be
        // reimplemented.
        // --------------------------------------------------------------------
        // if (mode_cols) {
        //     float temp_ref = vec_val[i];
        //     bool iscolumn = true;
        //     for (int j = 0; j != n; ++j) {
        //         if (temp_vec[j] != temp_ref) {
        //             iscolumn = false;
        //             break;
        //         }
        //     }
        //
        //     if (iscolumn && temp_ref > 0) {
        //         for (int j = 0; j != n; ++j) {
        //             *(nii_output_data + temp_vec_id[j]) = temp_ref;
        //             *(nii_output_extra_data + temp_vec_id[j]) = n;
        //         }
        //     }
        // }

        // --------------------------------------------------------------------
        // Find peak depth of maximum (Works with binary mask)
        // NOTE: temp_mask = 1 or 0
        // --------------------------------------------------------------------
        if (mode_peak) {
            float temp_ref = vec_val[i];
            float temp_max = vec_val[i];
            float temp_peak = vec_d[i];

            for (int j = 0; j != n; ++j) {
                if (temp_vec[j] > temp_max) {
                    temp_max = temp_vec[j];
                    temp_peak = temp_vec_d[j];
                }
            }
            *(nii_output_data + vec_voi_id[i]) = temp_peak;
            *(temp_nii_output_extra_data + vec_voi_id[i]) = static_cast<float>(n);
        }

        // --------------------------------------------------------------------
        // Count number of unique labels within each window
        // --------------------------------------------------------------------
        if (mode_count_uniques) {
            unordered_set<int> unique_elements;
            for (int j = 0; j != n; ++j) {
                unique_elements.insert(temp_vec[j]);
            }
            *(nii_output_data + vec_voi_id[i]) = unique_elements.size();
            *(temp_nii_output_extra_data + vec_voi_id[i]) = static_cast<float>(n);
        }
    }
    cout << endl;

    if (mode_median) {
        save_output_nifti(fout, "UVD_median_filter", nii_output, true);
    } else if (mode_min) {
        save_output_nifti(fout, "UVD_minpeaks", nii_output, true);
    } else if (mode_max) {
        save_output_nifti(fout, "UVD_max_filter", nii_output, true);
        save_output_nifti(fout, "UVD_max_filter_window_count", temp_nii_output_extra, true);
    } else if (mode_peak) {
        save_output_nifti(fout, "UVD_peak_filter", nii_output, true);
        save_output_nifti(fout, "UVD_peak_filter_window_count", temp_nii_output_extra, true);

    } else if (mode_cols) {
        save_output_nifti(fout, "UVD_columns_mode_filter", nii_output, true);
        save_output_nifti(fout, "UVD_columns_mode_filter_window_count_ratio", nii_output_extra, true);
        save_output_nifti(fout, "UVD_columns_mode_filter_window_count", temp_nii_output_extra, true);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
