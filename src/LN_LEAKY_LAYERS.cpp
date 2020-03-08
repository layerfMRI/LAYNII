

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_LEAKY_LAYERS: Layering algorithm based on iterative smothing.\n"
    "\n"
    "    This program derives 20 layers from GM/WM and GM/CSF border lines.\n"
    "\n"
    "Usage:\n"
    "    LN_LEAKY_LAYERS -rim rim.nii -dim 2 \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help.\n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -rim border         : Specify input dataset INT16.\n"
    "    -dim value          : Specify value (2 or 3) layer algorithm. \n"
    "                          Default is 3 (3D).\n"
    "\n"
    "Notes:\n"
    "    - This program now supports INT16, INT32 and FLOAT32. For other \n"
    "      data types use: \n"
    "          3dcalc -a rim.nii -datum short -expr 'a' -prefix rim_short.nii\n"
    "    - This can be 3D. Hence the rim file should be dmnii_smooth in all three\n"
    "      dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL, *fout = NULL;
    int ac, disp_float_eg = 0, dim;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options.
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-disp_float_example")) {
            disp_float_eg = 1;
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-dim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            dim = atof(argv[ac]);  // Assign pointer, no string copy
        } else {
            fprintf(stderr, " * * invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, " * * missing option '-rim'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_inputr = nifti_image_read(fin, 1);
    if (!nim_inputr) {
        fprintf(stderr, " * * failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_LEAKY_LAYERS");
    log_nifti_descriptives(nim_inputr);

    // Get dimensions of input
    int size_x = nim_inputr->nx;
    int size_y = nim_inputr->ny;
    int size_z = nim_inputr->nz;
    int size_t = nim_inputr->nt;
    int nr_voxels = nim_inputr->nvox;
    int nx = nim_inputr->nx;
    int nxy = nim_inputr->nx * nim_inputr->ny;
    int nxyz = nim_inputr->nx * nim_inputr->ny * nim_inputr->nz;
    float dX = nim_inputr->pixdim[1];
    float dY = nim_inputr->pixdim[2];
    float dZ = nim_inputr->pixdim[3];
    if (dim == 2) {
        dZ = 1000.;
        cout << "  Calculating layers only in 2D." << endl;
    }

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nim_input = copy_nifti_as_float32(nim_inputr);
    float* nim_input_data = static_cast<float*>(nim_input->data);

    // Allocate new nifti images
    nifti_image* nii_smooth = copy_nifti_as_float32(nim_input);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);
    nifti_image* nii_gaussw = copy_nifti_as_float32(nim_input);
    float* nii_gaussw_data = static_cast<float*>(nii_gaussw->data);
    nifti_image* nii_layer = copy_nifti_as_float32(nim_input);
    float* nii_layer_data = static_cast<float*>(nii_layer->data);
    nifti_image* nii_leaky = copy_nifti_as_float32(nim_input);
    float* nii_leaky_data = static_cast<float*>(nii_leaky->data);

    // ========================================================================
    cout << "  Debug 2 " << endl;
    float kernel_size = 2;  // Corresponds to one voxel sice.
    int vic = max(1., 2. * kernel_size);  // Ignore If voxel is too far
    float dist_i = 0.;
    cout << "  vic " << vic << endl;
    cout << "  kernel_size " << kernel_size << endl;

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_input_data + i) == 1) {
            *(nii_layer_data + i) = -200.;
        }
        if (*(nim_input_data + i) == 2) {
            *(nii_layer_data + i) = 200.;
        }
        if (*(nim_input_data + i) == 3) {
            *(nii_layer_data + i) = 0.;
        }
    }

    ///////////////////////////////
    // Start iterative loop here //
    ///////////////////////////////
    int iter_max = 400;

    for (int iter = 0; iter < iter_max; ++iter) {
        cout << "\r  Iteration: " << iter << " of " << iter_max << flush;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y-0; ++ix) {
                    int voxel_i = nxy * iz + nx * ix + iy;
                    *(nii_gaussw_data + voxel_i) = 0;
                    *(nii_smooth_data + voxel_i) = 0;
                    if (*(nim_input_data + voxel_i)  > 0) {
                        for (int iz_i = max(0, iz - vic);
                             iz_i < min(iz + vic + 1, size_y); ++iz_i) {
                            for (int iy_i = max(0, iy - vic);
                                 iy_i < min(iy + vic + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, ix - vic);
                                     ix_i < min(ix + vic + 1, size_y); ++ix_i) {
                                    int voxel_j = nxy * iz_i + nx * ix_i + iy_i;
                                    if (*(nim_input_data + voxel_j)  > 0) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz,
                                                      (float)ix_i, (float)iy_i, (float)iz_i,
                                                      dX, dY, dZ);
                                        if (dist_i < vic
                                            && *(nim_input_data + voxel_i) == 3) {
                                            *(nii_smooth_data + voxel_i) += *(nii_layer_data + voxel_j) * gaus(dist_i, kernel_size);
                                            *(nii_gaussw_data + voxel_i) += gaus(dist_i, kernel_size);
                                        }
                                    }
                                }
                            }
                        }
                        if (*(nii_gaussw_data + voxel_i) >= 0) {
                            *(nii_smooth_data + voxel_i) /= *(nii_gaussw_data + voxel_i);
                        }
                    }
                }
            }
        }
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    int voxel_i = nxy * iz + nx * ix + iy;
                    if (*(nim_input_data + voxel_i) == 1) {
                        *(nii_layer_data + voxel_i) = -200;
                    }
                    if (*(nim_input_data + voxel_i) == 2) {
                        *(nii_layer_data + voxel_i) = 200;
                    }
                    if (*(nim_input_data + voxel_i) == 3) {
                        *(nii_layer_data + voxel_i) = *(nii_smooth_data + voxel_i);
                    }
                }
            }
        }
    }  // Iteration loop closed


    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_leaky_data + i) = 2 + 19 * (*(nii_smooth_data + i) - (-200))
                                / (200 - (-200));
        if (*(nim_input_data + i) == 1) {
            *(nii_leaky_data + i) = 1;
        }
        if (*(nim_input_data + i) == 2) {
            *(nii_leaky_data + i) = 21;
        }
        if (*(nim_input_data + i) == 0) {
            *(nii_leaky_data + i) = 0;
        }
        if (*(nii_leaky_data + i) < 0) {
            *(nii_leaky_data + i) = 0;
        }
    }

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_input_data + i) != 0) {
            *(nii_leaky_data + i) = floor(22 - *(nii_leaky_data + i));
        }
        if (*(nii_leaky_data + i) == 20) {
            *(nii_leaky_data + i) = 19;
        }
        if (*(nii_leaky_data + i) == 21) {
            *(nii_leaky_data + i) = 20;
        }
    }
    save_output_nifti(fin, "leaky_layers", nii_leaky, true);

    cout << "  Finished." << endl;
    return 0;
}
