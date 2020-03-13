

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
    "    -help   : Show this help.\n"
    "    -rim    : Specify input dataset.\n"
    "    -dim    : Specify value (2 or 3) layer algorithm.Default is 3 (3D).\n"
    "\n"
    "Notes:\n"
    "    - This can be 3D. Hence the rim file should be dmsmooth in all\n"
    "      three dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac, dim;
    if (argc < 2) return show_help();

    // Process user options.
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-dim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            dim = atof(argv[ac]);
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
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, " * * failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_LEAKY_LAYERS");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_t = nii_input->nt;
    int nr_voxels = nii_input->nvox;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;
    float dX = nii_input->pixdim[1];
    float dY = nii_input->pixdim[2];
    float dZ = nii_input->pixdim[3];
    if (dim == 2) {
        dZ = 1000.;
        cout << "  Calculating layers only in 2D." << endl;
    }

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii_rim = copy_nifti_as_float32(nii_input);
    float* nii_rim_data = static_cast<float*>(nii_rim->data);

    // Allocate new nifti images
    nifti_image* smooth = copy_nifti_as_float32(nii_rim);
    float* smooth_data = static_cast<float*>(smooth->data);
    nifti_image* layers = copy_nifti_as_float32(nii_rim);
    float* layers_data = static_cast<float*>(layers->data);
    nifti_image* leaky = copy_nifti_as_float32(nii_rim);
    float* leaky_data = static_cast<float*>(leaky->data);

    // ========================================================================
    float kernel_size = 2;  // Corresponds to one voxel sice.
    int vic = max(1., 2. * kernel_size);  // Ignore if voxel is too far
    float d = 0.;
    cout << "  vic = " << vic << endl;
    cout << "  Kernel size = " << kernel_size << endl;

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 1) {
            *(layers_data + i) = -200.;
        }
        if (*(nii_rim_data + i) == 2) {
            *(layers_data + i) = 200.;
        }
        if (*(nii_rim_data + i) == 3) {
            *(layers_data + i) = 0.;
        }
    }

    // ========================================================================
    // Start iterative loop here
    // ========================================================================
    int iter_max = 400;

    for (int iter = 0; iter < iter_max; ++iter) {
        cout << "\r  Iteration: " << iter << " of " << iter_max << flush;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;
                    *(smooth_data + voxel_i) = 0;

                    if (*(nii_rim_data + voxel_i)  > 0) {
                        float total_weigth = 0;

                        int jz_start = max(0, iz - vic);
                        int jz_stop = min(iz + vic + 1, size_z);
                        int jy_start = max(0, iy - vic);
                        int jy_stop = min(iy + vic + 1, size_y);
                        int jx_start = max(0, ix - vic);
                        int jx_stop = min(ix + vic + 1, size_x);

                        for (int jz = jz_start; jz < jz_stop; ++jz) {
                            for (int jy = jy_start; jy < jy_stop; ++jy) {
                                for (int jx = jx_start; jx < jx_stop; ++jx) {
                                    int voxel_j = nxy * jz + nx * jy + jx;
                                    if (*(nii_rim_data + voxel_j)  > 0) {
                                        d = dist((float)ix, (float)iy, (float)iz,
                                                 (float)jx, (float)jy, (float)jz,
                                                 dX, dY, dZ);
                                        if (d < vic && *(nii_rim_data + voxel_i) == 3) {
                                            float w = gaus(d, kernel_size);
                                            *(smooth_data + voxel_i) += *(layers_data + voxel_j) * w;
                                            total_weigth += w;
                                        }
                                    }
                                }
                            }
                        }
                        if (total_weigth >= 0) {
                            *(smooth_data + voxel_i) /= total_weigth;
                        }
                    }
                }
            }
        }
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;
                    if (*(nii_rim_data + voxel_i) == 1) {
                        *(layers_data + voxel_i) = -200;
                    }
                    if (*(nii_rim_data + voxel_i) == 2) {
                        *(layers_data + voxel_i) = 200;
                    }
                    if (*(nii_rim_data + voxel_i) == 3) {
                        *(layers_data + voxel_i) = *(smooth_data + voxel_i);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    for (int i = 0; i < nr_voxels; ++i) {
        *(leaky_data + i) = 2 + 19 * (*(smooth_data + i) - (-200))
                                / (200 - (-200));
        if (*(nii_rim_data + i) == 1) {
            *(leaky_data + i) = 1;
        }
        if (*(nii_rim_data + i) == 2) {
            *(leaky_data + i) = 21;
        }
        if (*(nii_rim_data + i) == 0) {
            *(leaky_data + i) = 0;
        }
        if (*(leaky_data + i) < 0) {
            *(leaky_data + i) = 0;
        }
    }

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_rim_data + i) != 0) {
            *(leaky_data + i) = floor(22 - *(leaky_data + i));
        }
        if (*(leaky_data + i) == 20) {
            *(leaky_data + i) = 19;
        }
        if (*(leaky_data + i) == 21) {
            *(leaky_data + i) = 20;
        }
    }
    save_output_nifti(fin, "leaky_layers", leaky, true);

    cout << "  Finished." << endl;
    return 0;
}
