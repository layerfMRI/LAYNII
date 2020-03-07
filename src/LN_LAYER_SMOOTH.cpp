

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_LAYER_SMOOTH : Layering algorithm based on iterative smothing.\n"
    "\n"
    "    This program smooths data within layer or columns. In order to \n"
    "    avoid smoothing across masks, a crawler smooths only across \n"
    "    connected voxels.\n"
    "\n"
    "Usage:\n"
    "    LN_LAYER_SMOOTH -layer_file layers.nii -input activity_map.nii -FWHM 1\n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -layer_file : Nifti (.nii) file that contains layer or column masks.\n"
    "    -input      : Nifti (.nii) file that should be nii_smooth. It \n"
    "                  should have same dimensions as layer file.\n"
    "    -FWHM       : The amount of smoothing in mm.\n"
    "    -twodim     : (Optional) Smooth only in 2D. \n"
    "    -mask       : (Optional) Mask activity outside of layers. \n"
    "    -sulctouch  : (Optional) Allows smoothing across sucli. This is \n"
    "                  necessary, when you do heavy smoothing well bevond \n"
    "                  the spatial scale of the cortical thickness, or heavy\n"
    "                  curvature. It will make things slower. Note that this \n"
    "                  is best done with not too many layers. Otherwise a \n"
    "                  single layer has holes and is not connected.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* f_input = NULL, *f_layer = NULL;
    int ac, twodim = 0, do_masking = 0, sulctouch = 0;
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
            f_layer = argv[ac];  // No string copy, just pointer assignment
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
            f_input = argv[ac];  // No string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << "Smooth only in 2D."  << endl;
        } else if (!strcmp(argv[ac], "-sulctouch")) {
            sulctouch = 1;
            cout << "Smooth across sluci, might take longer."  << endl;
        } else if (!strcmp(argv[ac], "-mask")) {
            do_masking = 1;
            cout << "Set voxels to zero outside layers (mask option)"  << endl;
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

    // Read inputs including data
    nifti_image* nii_input = nifti_image_read(f_input, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_input);
        return 2;
    }

    nifti_image* nii_layer = nifti_image_read(f_layer, 1);
    if (!nii_layer) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_layer);
        return 2;
    }

    log_welcome("LN_LAYER_SMOOTH");
    log_nifti_descriptives(nii_input);
    log_nifti_descriptives(nii_layer);

    // Get dimensions of input
    const int size_z = nii_layer->nz;
    const int size_x = nii_layer->nx;
    const int size_y = nii_layer->ny;
    const int size_t = nii_layer->nt;
    const int nx = nii_layer->nx;
    const int nxy = nii_layer->nx * nii_layer->ny;
    const int nxyz = nii_layer->nx * nii_layer->ny * nii_layer->nz;
    const int nr_voxels = size_t * size_z * size_y * size_x;
    const float dX = nii_layer->pixdim[1];
    const float dY = nii_layer->pixdim[2];
    float dZ = nii_layer->pixdim[3];
    if (twodim == 1) dZ = 1000 * nii_layer->pixdim[3];


    // ========================================================================
    // Fix datatype issues
    nii_input = copy_nifti_as_float32(nii_input);
    float *nii_input_data = static_cast<float*>(nii_input->data);
    nii_input = copy_nifti_as_int32(nii_input);
    int32_t *nii_layer_data = static_cast<int32_t*>(nii_layer->data);

    // Allocate new niftis
    nifti_image* nii_smooth = copy_nifti_as_float32(nii_input);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);
    nifti_image* nii_gaussw = copy_nifti_as_float32(nii_input);
    float* nii_gaussw_data = static_cast<float*>(nii_gaussw->data);
    // ========================================================================

    // float kernel_size = 10;  // Corresponds to one voxel slice
    int vinc = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
    float dist_i = 0.;
    cout << "  Vicinity = " <<  vinc <<  endl;
    cout << "  FWHM = " <<  FWHM_val <<  endl;

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
    if (sulctouch == 0) {
        cout << "  Smoothing in layer, not considering sulci." << endl;

        for (int nr_layers_i = 1; nr_layers_i <= nr_layers; ++nr_layers_i) {
            cout << "\r    " << nr_layers_i << " of " << nr_layers << flush;

            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_y; ++iy) {
                    for (int ix = 0; ix < size_x - 0; ++ix) {
                        *(nii_gaussw_data + nxy * iz + nx * iy + ix) = 0;

                        if (*(nii_layer_data + nxy * iz + nx * iy + ix) == nr_layers_i) {

                            for (int iz_i = max(0, iz-vinc); iz_i < min(iz + vinc + 1, size_z-1); ++iz_i) {
                                for (int iy_i = max(0, iy-vinc); iy_i < min(iy + vinc + 1, size_y - 1); ++iy_i) {
                                    for (int ix_i = max(0, ix-vinc); ix_i < min(ix + vinc + 1, size_x - 1); ++ix_i) {
                                        if (*(nii_layer_data + nxy * iz_i + nx * iy_i + ix_i) == nr_layers_i) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            // cout << "  Debug 4 " << gaus(dist_i ,FWHM_val) <<   endl;
                                            // cout << "  Debug 5 " << dist_i  <<   endl;
                                            // if (*(nim_input_data + nxy * iz + nx * iy + ix) == 3) cout << "debug  4b " << endl;
                                            // dummy = *(layer_data + nxy*iz_i + nx*ix_i + iy_i);
                                            *(nii_smooth_data + nxy * iz + nx * iy + ix) = *(nii_smooth_data + nxy * iz + nx * iy + ix) + *(nii_input_data + nxy * iz_i + nx * iy_i + ix_i) * gaus(dist_i, FWHM_val);
                                            *(nii_gaussw_data + nxy * iz + nx * iy + ix) = *(nii_gaussw_data + nxy * iz + nx * iy + ix) + gaus(dist_i, FWHM_val);
                                        }
                                    }
                                }
                            }
                            if (*(nii_gaussw_data + nxy * iz + nx * iy + ix) > 0) {
                                *(nii_smooth_data + nxy * iz + nx * iy + ix) =
                                    *(nii_smooth_data + nxy * iz + nx * iy + ix)
                                    / *(nii_gaussw_data + nxy * iz + nx * iy + ix);
                            }
                        }
                        if (*(nii_layer_data + nxy * iz + nx * iy + ix) <= 0) {
                            *(nii_smooth_data + nxy * iz + nx * iy + ix) =
                                *(nii_input_data + nxy * iz + nx * iy + ix);
                        }
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
        // Allocating local connected vincinity file
        nifti_image* hairy_brain = copy_nifti_as_int32(nii_layer);
        int* hairy_brain_data = static_cast<int*>(hairy_brain->data);
        hairy_brain->scl_slope = 1.;
        int vinc_steps = 1;

        // Making sure voxels are zero
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_smooth_data + i) = 0;
        }

        vinc = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
        int nr_layers_i = 0;  // Running index
        cout << "  vinc " << vinc << endl;
        cout << "  FWHM_val " << FWHM_val << endl;
        cout << "  Starting within sulcal smoothing..." <<  endl;

        // For estimation of time
        int nr_vox_to_go_across = 0, running_index = 0, pref_ratio = 0;
        for (int i = 0; i < nr_voxels; ++i) {
            if (*(nii_layer_data + i) > 1) {
                nr_vox_to_go_across++;
            }
        }

        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x - 0; ++ix) {
                    if (*(nii_layer_data + nxy * iz + nx * iy + ix) > 0) {
                        running_index++;
                        if ((running_index * 100) / nr_vox_to_go_across != pref_ratio) {
                            cout << "\r " << (running_index * 100) / nr_vox_to_go_across <<  "% " << flush;
                            pref_ratio = (running_index * 100) / nr_vox_to_go_across;
                        }
                        nr_layers_i = *(nii_layer_data + nxy * iz + nx * iy + ix);
                        *(nii_gaussw_data + nxy * iz + nx * iy + ix) = 0;

                        /////////////////////////////////////////////////
                        // Find area that is not from the other sulcus //
                        /////////////////////////////////////////////////
                        for (int iz_i = max(0, iz - vinc - vinc_steps); iz_i <= min(iz + vinc + vinc_steps, size_z - 1); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc - vinc_steps); iy_i <= min(iy + vinc + vinc_steps, size_y - 1); ++iy_i) {
                                for (int ix_i = max(0, ix - vinc - vinc_steps); ix_i <= min(ix + vinc + vinc_steps, size_x - 1); ++ix_i) {
                                    // cout << iz_i << " " << iy_i << "  " << ix_i << "  " <<  size_z-1 << " " << size_x-1 << "  " << size_x-1 << "  " << endl;
                                    *(hairy_brain_data + nxy * iz_i + nx * iy_i + ix_i) = 0;
                                }
                            }
                        }
                        *(hairy_brain_data + nxy*iz + nx*ix + iy) = 1;

                        // Growing into neigbouring voxels.
                        for (int K_= 0; K_< vinc; K_++) {
                            for (int iz_ii = max(0, iz - vinc); iz_ii <= min(iz + vinc, size_z - 1); ++iz_ii) {
                                for (int iy_ii = max(0,iy - vinc); iy_ii <= min(iy + vinc, size_x - 1); ++iy_ii) {
                                    for (int ix_ii = max(0, ix-vinc); ix_ii <= min(ix + vinc, size_y - 1); ++ix_ii) {
                                        if (*(hairy_brain_data + nxy * iz_ii + nx * iy_ii + ix_ii) == 1) {
                                            for (int iz_i = max(0, iz_ii - vinc_steps); iz_i <= min(iz_ii + vinc_steps, size_z - 1); ++iz_i) {
                                                for (int iy_i = max(0, iy_ii - vinc_steps); iy_i<= min(iy_ii + vinc_steps, size_y - 1); ++iy_i) {
                                                    for (int ix_i = max(0, ix_ii - vinc_steps); ix_i <= min(ix_ii + vinc_steps, size_x - 1); ++ix_i) {
                                                        if (dist((float)ix_ii, (float)iy_ii, (float)iz_ii, (float)ix_i, (float)iy_i, (float)iz_i, 1, 1, 1) <= 1.74 && *(nii_layer_data + nxy * iz_i + nx * iy_i + ix_i) == nr_layers_i) {
                                                            *(hairy_brain_data + nxy * iz_i + nx * iy_i + ix_i) = 1;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Apply smoothing within each layer and within local patch
                        for (int iz_i = max(0, iz-vinc); iz_i <= min(iz + vinc, size_z - 1); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc); iy_i <= min(iy + vinc, size_y - 1); ++iy_i) {
                                for (int ix_i = max(0, ix-vinc); ix_i <= min(ix + vinc, size_x - 1); ++ix_i) {
                                    if (*(hairy_brain_data + nxy * iz_i + nx * iy_i + ix_i) == 1) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        *(nii_smooth_data + nxy * iz + nx * iy + ix) = *(nii_smooth_data + nxy * iz + nx * iy + ix) + *(nii_input_data + nxy * iz_i + nx * iy_i + ix_i) * gaus(dist_i, FWHM_val);
                                        *(nii_gaussw_data + nxy * iz + nx * iy + ix) = *(nii_gaussw_data + nxy * iz + nx * iy + ix) + gaus(dist_i, FWHM_val);
                                    }
                                }
                            }
                        }
                        if (*(nii_gaussw_data + nxy * iz + nx * iy + ix) > 0) *(nii_smooth_data + nxy * iz + nx * iy + ix) = *(nii_smooth_data + nxy * iz + nx * iy + ix) / *(nii_gaussw_data + nxy * iz + nx * iy + ix);
                    } // if scope  if (*(nii_layer_data + nxy * iz + nx * iy + ix) > 0){ closed

                }
            }
        }
        // for(int iz=0; iz<size_z; ++iz){
        //     for(int iy=0; iy<size_x; ++iy){
        //         for(int ix=0; ix<size_y; ++ix){
        //             if(*(nii_input_data + nxy * iz + nx * iy + ix) > 0) {
        //                 *(nii_input_data + nxy * iz + nx * iy + ix) = *(nii_smooth_data + nxy * iz + nx * iy + ix);
        //             }
        //         }
        //     }
        // }
        const char* fout_3 = "hairy_brain.nii";
        if (nifti_set_filenames(hairy_brain, fout_3, 1, 1)) return 1;
        nifti_image_write(hairy_brain);
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

    nii_smooth->scl_slope = nii_input->scl_slope;

    if (nii_input->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    // Output file name
    // const char* fout_4 = "leaky_layers.nii" ;
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) return 1;
    // nifti_image_write(leak_layer);

    // const char* fout_5 = "input_file.nii" ;
    // if (nifti_set_filenames(nii_input, fout_5, 1, 1)) return 1;
    // nifti_image_write(nii_input);

    // const char* fout_2 = "mask.nii" ;
    // if (nifti_set_filenames(nii_layer, fout_2, 1, 1)) return 1;
    // nifti_image_write(nii_layer);

    string prefix = "nii_smooth_";
    string filename = (string) (f_input);
    string outfilename = prefix + filename;
    log_output(outfilename.c_str());

    const char* fout_1 = outfilename.c_str();
    if (nifti_set_filenames(nii_smooth, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(nii_smooth);

    cout << "  Finished." << endl;
    return 0;
}
