
// TODO(Faruk): Seems there might be an issue with the gaussian kernel's
// symmetry and size corresponding to what is written in CLI

// To do make the vincinity direction specific vinc_x, vinc_y, vinc_z


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_DEVEIN : Program that tries to remove the macrovascular component in BOLD.\n"
    "\n"
    "    This program tries to estimate the micrcovascular component in layer-fMRI GE-BOLD \n"
    "    It does so by using an estimate of the local macrovascular blood volume ALF \n"
    "\n"
    "\n"
    "Usage:\n"
    "    LN2_DEVEIN -layer_file layers.nii -column_file columns.nii -input input time_series.nii -ALF ALF.nii  -FWHM 1\n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -layer_file    : Nifti (.nii) file that contains layers. \n"
    "    -column_file   : Nifti (.nii) file that should be smooth. It \n"
    "    -FWHM          : optional amount of smoothing in mm.\n"
    "    -input         : BOLD file that should be corrected from macro-vascular contaminations\n"
    "                     This can be a time series or an activity map (not z-scores though). \n"
    "    -ALF           : File with estiates of Amplitude of low frequencies as a correlate ot venous CBV  \n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* f_input = NULL, *f_layer = NULL, *f_column = NULL, *f_ALF = NULL;
    int ac, do_masking = 0, sulctouch = 0;
    float FWHM_val = 0;
    if (argc < 3) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layer_file\n");
                return 1;
            }
            f_layer = argv[ac];
        } else if (!strcmp(argv[ac], "-column_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -column_file\n");
                return 1;
            }
            f_column = argv[ac];
        } else if (!strcmp(argv[ac], "-ALF")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -ALF\n");
                return 1;
            }
            f_ALF = argv[ac];
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);  // No string copy, pointer assignment
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            f_input = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!f_input) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!f_layer) {
        fprintf(stderr, "** missing option '-layer_file'\n");
        return 1;
    }
    if (!f_column) {
        fprintf(stderr, "** missing option '-column_file'\n");
        return 1;
    }
    if (!f_ALF) {
        fprintf(stderr, "** missing option '-ALF'\n");
        return 1;
    }

    // Read inputs including data
    nifti_image* nii = nifti_image_read(f_input, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_input);
        return 2;
    }

    nifti_image* nii_layeri = nifti_image_read(f_layer, 1);
    if (!nii_layeri) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_layer);
        return 2;
    }
    
    nifti_image* nii_columni = nifti_image_read(f_column, 1);
    if (!nii_columni) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_column);
        return 2;
    }
    
    nifti_image* nii_ALFi = nifti_image_read(f_ALF, 1);
    if (!nii_ALFi) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_ALF);
        return 2;
    }

    log_welcome("LN2_DEVEIN");
    log_nifti_descriptives(nii);
    log_nifti_descriptives(nii_layeri);
    log_nifti_descriptives(nii_columni);
    log_nifti_descriptives(nii_ALFi);

    // Get dimensions of input
    const int size_z = nii->nz;
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int nx = nii->nx;
    const int nxy = nii->nx * nii->ny;
    const int nr_voxels = size_z * size_y * size_x;
    const float dX = nii->pixdim[1];
    const float dY = nii->pixdim[2];
    const float dZ = nii->pixdim[3];

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii);
    float *nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* nii_layer = copy_nifti_as_int32(nii_layeri);
    int32_t *nii_layer_data = static_cast<int32_t*>(nii_layer->data);
    nifti_image* nii_column = copy_nifti_as_int32(nii_columni);
    int32_t *nii_column_data = static_cast<int32_t*>(nii_column->data);
    nifti_image* nii_ALF = copy_nifti_as_float32(nii_ALFi);
    float *nii_ALF_data = static_cast<float*>(nii_ALF->data);

    // Allocate new niftis
    nifti_image *nii_smooth = copy_nifti_as_float32(nii_input);
    float *nii_smooth_data = static_cast<float*>(nii_smooth->data);
    nifti_image* nii_gaussw = copy_nifti_as_float32(nii_input);
    float *nii_gaussw_data = static_cast<float*>(nii_gaussw->data);
/*
    // Zero new images
    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_smooth_data + i) = 0;
        *(nii_gaussw_data + i) = 0;
    }

    // ========================================================================
    // TODO(Faruk): Why using dX but not others? Need to ask Renzo about this.
    int vic = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
    cout << "  Vicinity = " << vic << endl;
    cout << "  FWHM = " << FWHM_val << endl;

    ///////////////////////////
    // Find number of layers //
    ///////////////////////////
    int32_t nr_layers = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_layer_data + i) > nr_layers) {
            nr_layers = *(nii_layer_data + i);
        }
    }
    cout << "  There are " << nr_layers << " layers to smooth within." << endl;

    ////////////////////
    // SMOOTHING LOOP //
    ////////////////////
    // For time estimation
    int nr_vox_to_loop = 0, idx = 0, prev_n = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_layer_data + i) > 0) {
            nr_vox_to_loop++;
        }
    }

    if (sulctouch == 0) {
        cout << "  Smoothing in layer, not considering sulci." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;
                    int layer_i = *(nii_layer_data + voxel_i);

                    if (layer_i > 0) {
                        idx += 1;
                        int n = (idx * 100) / nr_vox_to_loop;
                        if (n != prev_n) {
                            cout << "\r    " << n <<  "%" << flush;
                            prev_n = n;
                        }

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
                                    if (*(nii_layer_data + voxel_j) == layer_i) {
                                        float d = dist((float)ix, (float)iy, (float)iz,
                                                       (float)jx, (float)jy, (float)jz,
                                                       dX, dY, dZ);
                                        float g = gaus(d, FWHM_val);
                                        *(nii_smooth_data + voxel_i) += *(nii_input_data + voxel_j) * g;
                                        *(nii_gaussw_data + voxel_i) += g;
                                    }
                                }
                            }
                        }
                        // Normalize
                        *(nii_smooth_data + voxel_i) /= *(nii_gaussw_data + voxel_i);
                    } else {
                        *(nii_smooth_data + voxel_i) = *(nii_input_data + voxel_i);
                    }
                }
            }
        }
        cout << endl;
    }

    ///////////////////////////////////////////////////////
    // if requested, smooth only within connected layers //
    ///////////////////////////////////////////////////////
    if (sulctouch == 1) {
        // Allocating local connected vicinity file
        nifti_image* hairy_brain = copy_nifti_as_int32(nii_layer);
        int32_t* hairy_brain_data = static_cast<int32_t*>(hairy_brain->data);
        hairy_brain->scl_slope = 1.;
        int vic_steps = 1;

        vic = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
        cout << "  vic " << vic << endl;
        cout << "  FWHM_val " << FWHM_val << endl;
        cout << "  Starting within sulcus smoothing..." <<  endl;

        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;

                    if (*(nii_layer_data + voxel_i) > 0) {
                        idx++;
                        int n = (idx * 100) / nr_vox_to_loop;
                        if (n != prev_n) {
                            cout << "\r " << n <<  "% " << flush;
                            prev_n = n;
                        }
                        int layer_i = *(nii_layer_data + voxel_i);

                        /////////////////////////////////////////////////
                        // Find area that is not from the other sulcus //
                        /////////////////////////////////////////////////
                        int jz_start = max(0, iz - vic - vic_steps);
                        int jz_stop = min(iz + vic + vic_steps, size_z - 1);
                        int jy_start = max(0, iy - vic - vic_steps);
                        int jy_stop = min(iy + vic + vic_steps, size_y - 1);
                        int jx_start = max(0, ix - vic - vic_steps);
                        int jx_stop = min(ix + vic + vic_steps, size_x - 1);

                        for (int jz = jz_start; jz <= jz_stop; ++jz) {
                            for (int jy = jy_start; jy <= jy_stop; ++jy) {
                                for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                    *(hairy_brain_data + nxy * jz + nx * jy + jx) = 0;
                                }
                            }
                        }
                        *(hairy_brain_data + voxel_i) = 1;

                        // Grow into neigbouring voxels.
                        for (int K_= 0; K_< vic; K_++) {
                            int kz_start = max(0, iz - vic);
                            int kz_stop = min(iz + vic, size_z - 1);
                            int ky_start = max(0, iy - vic);
                            int ky_stop = min(iy + vic, size_x - 1);
                            int kx_start = max(0, ix - vic);
                            int kx_stop = min(ix + vic, size_y - 1);

                            for (int kz = kz_start; kz <= kz_stop; ++kz) {
                                for (int ky = ky_start; ky <= ky_stop; ++ky) {
                                    for (int kx = kx_start; kx <= kx_stop; ++kx) {
                                        if (*(hairy_brain_data + nxy * kz + nx * ky + kx) == 1) {
                                            int mz_start = max(0, kz - vic_steps);
                                            int mz_stop = min(kz + vic_steps, size_z - 1);
                                            int my_start = max(0, ky - vic_steps);
                                            int my_stop = min(ky + vic_steps, size_y - 1);
                                            int mx_start = max(0, kx - vic_steps);
                                            int mx_stop = min(kx + vic_steps, size_x - 1);

                                            for (int mz = mz_start; mz <= mz_stop; ++mz) {
                                                for (int my = my_start; my <= my_stop; ++my) {
                                                    for (int mx = mx_start; mx <= mx_stop; ++mx) {
                                                        if (dist((float)kx, (float)ky, (float)kz, (float)mx, (float)my, (float)mz, 1, 1, 1) <= 1.74
                                                            && *(nii_layer_data + nxy * mz + nx * my + mx) == layer_i) {
                                                            *(hairy_brain_data + nxy * mz + nx * my + mx) = 1;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Smooth within each layer and within local patch
                        jz_start = max(0, iz - vic);
                        jz_stop = min(iz + vic, size_z - 1);
                        jy_start = max(0, iy - vic);
                        jy_stop = min(iy + vic, size_y - 1);
                        jx_start = max(0, ix - vic);
                        jx_stop = min(ix + vic, size_x - 1);

                        for (int jz = jz_start; jz <= jz_stop; ++jz) {
                            for (int jy = jy_start; jy <= jz_stop; ++jy) {
                                for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                    if (*(hairy_brain_data + nxy * jz + nx * jy + jx) == 1) {
                                        float d = dist((float)ix, (float)iy, (float)iz,
                                                       (float)jx, (float)jy, (float)jz,
                                                       dX, dY, dZ);
                                        float g = gaus(d, FWHM_val);

                                        *(nii_smooth_data + voxel_i) += *(nii_input_data + nxy * jz + nx * jy + jx) * g;
                                        *(nii_gaussw_data + voxel_i) += g;
                                    }
                                }
                            }
                        }
                        if (*(nii_gaussw_data + voxel_i) > 0) {
                        *(nii_smooth_data + voxel_i) /= *(nii_gaussw_data + voxel_i);
                        }
                    }

                }
            }
        }
        save_output_nifti(f_input, "hairy_brain", hairy_brain, false);
    }
    cout << "  Smoothing is done. " <<  endl;

    ////////////////////////////////
    // Masking if it is it wanted //
    ////////////////////////////////
    if (do_masking == 1) {
        for (int i = 0; i < nr_voxels; ++i)
            if (*(nii_layer_data + i) == 0) {
                *(nii_smooth_data + i) = 0;
        }
    }
*/
    save_output_nifti(f_input, "deconvolved", nii_smooth, true);

    cout << "  Finished." << endl;
    return 0;
}
