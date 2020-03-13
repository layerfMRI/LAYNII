

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_GRADSMOOTH : Local gradient based smoothing.\n"
    "\n"
    "Usage:\n"
    "    LN_GRADSMOOTH -input activity_map.nii -gradfile gradfile.nii -FWHM 1 -within -selectivity 0.1 \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -input       : Nifti (.nii) that will be smoothed.\n"
    "    -gradfile    : Nifti (.nii) used to estimate local gradients.\n"
    "                   Only the first time point of this file is used. It \n"
    "                   should have the same spatial dimensions as the input.\n"
    "    -FWHM        : Amount of smoothing in mm.\n"
    "    -twodim      : (Optional) Smooth in 2 dimensions only. \n"
    "    -mask        : (Optional) Nifti (.nii) that is used mask activity \n"
    "                   outside. This option can speed up processing.\n"
    "    -within      : (Optional) Determines that smoothing should happen \n"
    "                   within similar values, not across different values.\n"
    "    -acros       : (Optional) Determines that smoothing should happen \n"
    "                   across different values, not within similar values.\n"
    "                   NOTE: This option is not working yet.\n"
    "    -selectivity : (Optional) Makes the smoothing more or less specific \n"
    "                   to a certain gradient range. 0.05 is only within \n"
    "                   very similar values. 0.9 is almost independent of \n"
    "                   the gradient file 0.1 is default.\n"
    "\n"
    "Notes:\n"
    "    - If you run this on EPI-T1 data consider making them pretty, E.g: \n"
    "        start_bias_field.sh T1.nii \n"
    "        denoise_me.sh bico_T1.nii \n"
    "        short_me.sh denoised_bico_T1.nii \n"
    "        smooth_me.sh denoised_bico_T1.nii 0.5 \n"
    "        mv smooth_denoised_bico_T1.nii new_T1.nii \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * fin_2 = NULL, * fin_1 = NULL, * fin_3 = NULL;
    int ac, twodim = 0, do_masking = 0, within = 0, acros = 0;
    float FWHM_val = 0, selectivity = 0.1;
    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-gradfile")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -gradfile\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << "Smooth only in 2D." << endl;
        } else if (!strcmp(argv[ac], "-within")) {
            within = 1;
            cout << "Within similar values." << endl;
        } else if (!strcmp(argv[ac], "-acros")) {
            acros = 1;
            cout << "Across different values." << endl;
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            fin_3 = argv[ac];
            do_masking = 1;
            cout << "Set voxels to 0 outside layers (mask option)." << endl;
        } else if (!strcmp(argv[ac], "-selectivity")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -selectivity\n");
                return 1;
            }
            selectivity = atof(argv[ac]);
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-gradfile'\n");
        return 1;
    }
    if (acros + within !=1) {
        cout << "Please choose within or across." << endl;
        return 2;
    }
    if (acros ==1) {
        cout << "Smoothing across gradients is not implemented yet. Use -within instead  " << endl;
        return 2;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_GRADSMOOTH");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    int size_x = nii2->nx;
    int size_y = nii2->ny;
    int size_z = nii2->nz;
    int size_t = nii1->nt;
    int nx = nii2->nx;
    int nxy = nii2->nx * nii2->ny;
    int nxyz = nii2->nx * nii2->ny * nii2->nz;
    float dX = nii2->pixdim[1];
    float dY = nii2->pixdim[2];
    float dZ = nii2->pixdim[3];
    // If you are running the smoothing in 2D, it will still go thought the
    // entire pipeline the only difference is that the weights in a certain
    // direction are suppressed doing it in 2D, will not speed up the program
    if  (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* nii_grad = copy_nifti_as_float32(nii2);
    float* nii_grad_data = static_cast<float*>(nii_grad->data);

    // Allocate new nifti images
    nifti_image* smooth = copy_nifti_as_float32(nii_input);
    float* smooth_data = static_cast<float*>(smooth->data);

    // Mask nifti related part
    nifti_image* nii_mask;
    float* nii_mask_data;
    if (do_masking == 1) {
        nifti_image* nii3 = nifti_image_read(fin_3, 1);
        if (!nii3) {
            fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_3);
            return 2;
        }
        log_nifti_descriptives(nii3);
        nii_mask = copy_nifti_as_int32(nii3);
        nii_mask_data = static_cast<float*>(nii_mask->data);
    } else {  // When no mask, make all voxels 1
        nii_mask = copy_nifti_as_int32(nii1);
        nii_mask_data = static_cast<float*>(nii_mask->data);
        for (int i = 0; i < nii2->nvox; ++i) {
            *(nii_mask_data + i) = 1;
        }
    }
    // ========================================================================

    cout << "  Time dimension of smooth. Output file: " << smooth->nt << endl;

    int vic = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
    cout << "  vic " << vic << endl;
    cout << "  FWHM_val " << FWHM_val << endl;

    // ========================================================================
    // Finding the range of gradient values
    // ========================================================================

    // Values that I need to characterize the local signals in the vicinity.
    int nr_vox_in_vic = (2 * vic + 1) * (2 * vic + 1) * (2 * vic + 1);
    double vec1[nr_vox_in_vic];
    for (int it = 0; it < nr_vox_in_vic; it++) {
        vec1[it] = 0;
    }

    // Estimate and output of program process and how much longer it will take.
    int nr_voxels = size_z * size_x * size_y;
    int running_index = 0, pref_ratio = 0;
    if (do_masking == 1) {
        nr_voxels = 0;
        for (int i = 0; i < nii2->nvox; ++i) {
            if (*(nii_mask_data + i) > 0) {
                nr_voxels = nr_voxels + 1;
            }
        }
    }
    cout << "  The number of voxels to go across = " << nr_voxels << endl;

    // ========================================================================
    // Smoothing loop
    // ========================================================================
    cout << "  Big smoothing loop is being done..." << endl;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                if (*(nii_mask_data + voxel_i) != 0) {
                    // Write out how many more voxels I have to go through.
                    running_index++;
                    if ((running_index * 100) / nr_voxels != pref_ratio) {
                        cout << "\r  Progress: "
                             << (running_index * 100) / nr_voxels
                             << "%" << flush;
                        pref_ratio = (running_index * 100)/nr_voxels;
                    }
                    // Clean kitchen
                    float w = 0;
                    *(smooth_data + voxel_i) = 0;
                    nr_vox_in_vic = 0;
                    float val_i = *(nii_grad_data + voxel_i);

                    // Examining the environment and determining what the
                    // signal intensities are and what its distribution are
                    int jz_start = max(0, iz - vic);
                    int jz_stop = min(iz + vic, size_z - 1);
                    int jy_start = max(0, iy - vic);
                    int jy_stop = min(iy + vic, size_y - 1);
                    int jx_start = max(0, ix - vic);
                    int jx_stop = min(ix + vic, size_x - 1);
                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;
                                float val_j = *(nii_grad_data + voxel_j);
                                vec1[nr_vox_in_vic] = static_cast<double>(val_j);
                                nr_vox_in_vic++;
                            }
                        }
                    }

                    // Standard deviation of the signal valued in the vicinity.
                    // This is necessary to normalize how many voxels are
                    // contributing to the local smoothing.
                    float grad_stdev = static_cast<float>(ren_stdev(vec1, nr_vox_in_vic));

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;
                                float val_j = *(nii_grad_data + voxel_j);

                                float d1 = dist((float)ix, (float)iy, (float)iz,
                                                (float)jx, (float)jy, (float)jz,
                                                dX, dY, dZ);
                                float g1 = gaus(d1, FWHM_val);

                                float d2 = abs(val_i - val_j);
                                float g2 = gaus(d2, grad_stdev * selectivity);

                                float g3 = gaus(0, grad_stdev * selectivity);

                                float temp_w =  g1 * g2 / g3;

                                // NOTE(Renzo): Gauss data are important to
                                // kernel size changes (e.g. at edge of images).
                                // avoid local scaling differences. When the
                                // This is a geometric parameter and only need
                                // to be calculated for one time point. This
                                // might be avoidable, if the Gauss fucnction
                                // is better normalized.
                                w += temp_w;
                                for (int it = 0; it < size_t; ++it) {
                                    *(smooth_data + nxyz * it + voxel_i) +=
                                        *(nii_input_data + nxyz * it + voxel_j) * temp_w;
                                }
                            }
                        }
                    }
                    // Scaling signal intensity with the overall Gauss leakage
                    if (w > 0) {
                        for (int it = 0; it < size_t; ++it) {
                            *(smooth_data + nxyz * it + voxel_i) /= w;
                        }
                    }
                }
            }
        }
    }
    cout << endl;
    save_output_nifti(fin_1, "smooth", smooth, true);

    cout << "  Finished." << endl;
    return 0;
}
