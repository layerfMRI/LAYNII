
#include "../dep/laynii_lib.h"
#include <limits>


int show_help(void) {
    printf(
    "LN2_RIM_POLISH: Smooth the cortical gray matter borders. Default parameters\n"
    "                are optimized for 0.2 mm isotropic images."
    "\n"
    "Usage:\n"
    "    LN2_RIM_POLISH -rim rim.nii.gz \n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -rim           : A segmented image. This image must use 1 to code outer\n"
    "                     gray matter border voxels (facing mostly CSF), 2 to code inner\n"
    "                     gray matter border voxels (facing mostly white matter), and\n"
    "                     3 to code pure gray matter voxels.\n"
    "    -steps         : (Optional) Number of erosion dilation steps applied before\n"
    "                     the Gaussian smoothing step. It is recommended to only\n"
    "                     increase this for <0.2 mm isotropic resolution \n"
    "                     corticalimages. '1' by default, chosen for 0.2 mm iso. images.\n"
    "    -iter_smooth   : (Optional) Number of smoothing iterations. Higher values\n"
    "                     will result in smoother boundaries. However, this also means\n"
    "                     that fundi of the sulci and crowns of the gyri might get \n"
    "                     smoothed out. '6' by default, chosen for 0.2 mm iso. images.\n"
    "    -steps_voronoi : (Optional) Number of voronoi dilation steps. Useful for\n"
    "                     preventing smoothing artifacts around the edges of partially\n" 
    "                     segmented volumes. '5' by default, chosen for 0.2 mm iso. images.\n"
    "    -debug         : (Optional) Save extra intermediate outputs.\n"
    "    -output        : (Optional) Output filename, including .nii or\n"
    "                     .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
    bool use_outpath = false, mode_debug=false;
    char *fin = NULL, *fout = NULL;
    int ac;
    int steps=1, iter_smooth=6, steps_voronoi=5;

    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -rim\n");
                return 1;
            }
            fin = argv[ac];
            fout = fin;
        } else if (!strcmp(argv[ac], "-steps")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -steps\n");
                return 2;
            }
            steps = std::stoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-iter_smooth")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -iter_smooth\n");
                return 2;
            }
            iter_smooth = std::stoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-steps_voronoi")) {
            steps_voronoi = std::stoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            use_outpath = true ;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }


    // Read input dataset
    nifti_image *nii_in = nifti_image_read(fin, 1);
    if (!nii_in) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN2_RIM_POLISH");
    log_nifti_descriptives(nii_in);

    // Get dimensions of input
    const int size_x = nii_in->nx;
    const int size_y = nii_in->ny;
    const int size_z = nii_in->nz;
    const int nr_voxels = size_z * size_y * size_x;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    // Prepare images
    nifti_image *nii_rim = copy_nifti_as_int16(nii_in);
    int16_t* nii_rim_data = static_cast<int16_t*>(nii_rim->data);

    // ========================================================================
    // Voronoi dilate initial segmentation (useful for partial segmentations)
    // ========================================================================
    nifti_image* nii_temp = copy_nifti_as_int16(nii_rim);
    int16_t* nii_temp_data = static_cast<int16_t*>(nii_temp->data);

    uint32_t ix, iy, iz, j;
    for (uint32_t n = 0; n != steps_voronoi; ++n) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_temp_data + i) == 0) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_rim_data + i) = *(nii_temp_data + j);
                    }
                }
            }
        }
        // Swap temp1 and temp2
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_temp_data + i) = *(nii_rim_data + i);
        }
    }

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_temp_data + i) = *(nii_rim_data + i);
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "voronoi_dilated", nii_rim, false);
    }

    // ========================================================================
    // Binarize wm and wmgm images
    // ========================================================================
    // Prepare temporary images
    nifti_image *nii_wm = copy_nifti_as_int16(nii_rim);
    int16_t* nii_wm_data = static_cast<int16_t*>(nii_wm->data);
    nifti_image *nii_wmgm = copy_nifti_as_int16(nii_rim);
    int16_t* nii_wmgm_data = static_cast<int16_t*>(nii_wmgm->data);

    cout << "  Binarizing white matter (wm) and white + gray matter (wmgm)..." << endl;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_wm_data + i) == 2) {
            *(nii_wm_data + i) = 1;
        } else {
            *(nii_wm_data + i) = 0;
        }
    }

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_wmgm_data + i) == 2 || *(nii_wmgm_data + i) == 3) {
            *(nii_wmgm_data + i) = 1;
        } else {
            *(nii_wmgm_data + i) = 0;
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wm", nii_wm, false);
        save_output_nifti(fout, "wmgm", nii_wmgm, false);
    }

    // ========================================================================
    // Dilate white matter
    // ========================================================================
    cout << "  Polishing white matter (wm)..." << endl;
    cout << "    Dilating..." << endl;

    // Clean temporary nifti
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_temp_data + i) = *(nii_wm_data + i);
    }

    // Dilate
    for (uint32_t n = 0; n != steps; ++n) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_temp_data + i) == 0) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
            }
        }
        // Swap temp1 and temp2
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_temp_data + i) = *(nii_wm_data + i);
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wm_dilated", nii_wm, false);
    }

    // ------------------------------------------------------------------------
    // Smooth dilated white matter
    // ------------------------------------------------------------------------
    cout << "    Smoothing..." << endl;

    nifti_image* nii_wm_smth = copy_nifti_as_float32(nii_wm);
    float* nii_wm_smth_data = static_cast<float*>(nii_wm_smth->data);

    // NOTE(Faruk): Mask is needed to constrain the iterative smoothing.
    nifti_image* nii_mask = copy_nifti_as_int16(nii_rim);
    int16_t* nii_mask_data = static_cast<int16_t*>(nii_mask->data);
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) != 0) {
            *(nii_mask_data + i) = 1;
        } else {
            *(nii_mask_data + i) = 0;
        }
    }

    nii_wm_smth = iterative_smoothing(nii_wm_smth, iter_smooth, nii_mask, 1);
    nii_wm_smth_data = static_cast<float*>(nii_wm_smth->data);

    if (mode_debug == true) {
        save_output_nifti(fout, "wm_dilated_smoothed", nii_wm_smth, false);
    }

    // ------------------------------------------------------------------------
    // Binarize smoothed wm
    // ------------------------------------------------------------------------
    cout << "    Binarizing..." << endl;

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_wm_smth_data + i) >= 0.5) {
            *(nii_wm_data + i) = 1;
        } else {
            *(nii_wm_data + i) = 0;
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wm_dilated_smoothed_binarized", nii_wm, false);
    }

    // ------------------------------------------------------------------------
    // Erode back wm
    // ------------------------------------------------------------------------
    cout << "    Eroding..." << endl;
    // Clean temporary nifti
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_temp_data + i) = *(nii_wm_data + i);
    }

    // Dilate
    for (uint32_t n = 0; n != steps; ++n) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_temp_data + i) == 1) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wm_data + i) = *(nii_temp_data + j);
                    }
                }
            }
        }
        // Swap temp1 and temp2
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_temp_data + i) = *(nii_wm_data + i);
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wm_dilated_smoothed_binarized_eroded", nii_wm, false);
    }

    // ========================================================================
    // Erode wmgm (white + gray matter)
    // ========================================================================
    cout << "  Polishing white + gray matter (wmgm)..." << endl;
    cout << "    Eroding..." << endl;

    // Clean temporary nifti
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_temp_data + i) = *(nii_wmgm_data + i);
    }

    // Dilate
    for (uint32_t n = 0; n != steps; ++n) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_temp_data + i) == 1) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_temp_data + j) == 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
            }
        }
        // Swap temp1 and temp2
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_temp_data + i) = *(nii_wmgm_data + i);
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wmgm_eroded", nii_wmgm, false);
    }

    // ------------------------------------------------------------------------
    // Smooth dilated wmgm
    // ------------------------------------------------------------------------
    cout << "    Smoothing..." << endl;

    nifti_image* nii_wmgm_smth = copy_nifti_as_float32(nii_wmgm);
    float* nii_wmgm_smth_data = static_cast<float*>(nii_wmgm_smth->data);

    nii_wmgm_smth = iterative_smoothing(nii_wmgm_smth, iter_smooth, nii_mask, 1);
    nii_wmgm_smth_data = static_cast<float*>(nii_wmgm_smth->data);

    if (mode_debug == true) {
        save_output_nifti(fout, "wmgm_eroded_smoothed", nii_wmgm_smth, false);
    }

    // ------------------------------------------------------------------------
    // Binarize smoothed wmgm
    // ------------------------------------------------------------------------
    cout << "    Binarizing..." << endl;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_wmgm_smth_data + i) >= 0.5) {
            *(nii_wmgm_data + i) = 1;
        } else {
            *(nii_wmgm_data + i) = 0;
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wmgm_eroded_smoothed_binarized", nii_wmgm, false);
    }

    // ------------------------------------------------------------------------
    // Dilate back wmgm
    // ------------------------------------------------------------------------
    cout << "    Eroding..." << endl;

    // Clean temporary nifti
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_temp_data + i) = *(nii_wmgm_data + i);
    }

    // Dilate
    for (uint32_t n = 0; n != steps; ++n) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_temp_data + i) == 0) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_temp_data + j) != 0) {
                        *(nii_wmgm_data + i) = *(nii_temp_data + j);
                    }
                }
            }
        }
        // Swap temp1 and temp2
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_temp_data + i) = *(nii_wmgm_data + i);
        }
    }

    if (mode_debug == true) {
        save_output_nifti(fout, "wmgm_eroded_smoothed_binarized_dilated", nii_wmgm, false);
    }

    // ========================================================================
    // Combine
    // ========================================================================
    cout << "  Combining..." << endl;

    nifti_image *nii_rim_orig = copy_nifti_as_int16(nii_in);
    int16_t* nii_rim_orig_data = static_cast<int16_t*>(nii_rim_orig->data);

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_orig_data + i) != 0) {
            if (*(nii_wm_data + i) == 1) {
                *(nii_temp_data + i) = 2;
            } else if (*(nii_wmgm_data + i) == 1) {
                *(nii_temp_data + i) = 3;
            } else if (*(nii_rim_orig_data + i) != 0) {
                *(nii_temp_data + i) = 1;
            } 
        } else {
                *(nii_temp_data + i) = 0;
        }
    }

    // Final output
    save_output_nifti(fout, "polished", nii_temp, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
