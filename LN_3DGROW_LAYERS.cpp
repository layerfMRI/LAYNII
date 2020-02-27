
#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_3DGROW_LAYERS: Cortical gray matter layering.\n"
    "\n"
    "Usage:\n"
    "    LN_3DGROW_LAYERS -rim rim.nii \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help. \n"
    "    -rim        : Specify input dataset.\n"
    "    -nr_layers  : Number of layers. Default is 3.\n"
    "\n"
    "Notes:\n"
    "    - Datatype of 'rim.nii' needs to be INT16.\n"
    "    - This is 3D. Hence rim.nii file should be dmsmooth in all 3D.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image*nii_rim = NULL;
    char* fin = NULL;
    int ac, nr_layers = 3;
    float column_size = 7;
    if (argc < 2) {
        return show_help();   // typing '-help' is sooo much work
    }

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
        } else if (!strcmp(argv[ac], "-nr_layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_layers\n");
            } else {
                nr_layers = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-column_size")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -column_size\n");
            } else {
                column_size = atof(argv[ac]);
            }
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }

    // Read input dataset, including data
    nii_rim = nifti_image_read(fin, 1);
    if (!nii_rim) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }
    int16_t *nii_rim_data = static_cast<__int16_t*>(nii_rim->data);

    log_welcome("LN_3DGROW_LAYERS");
    log_nifti_descriptives(nii_rim);

    cout << "  Nr. layers: " << nr_layers << endl;

    // Get dimensions of input
    const int size_z = nii_rim->nz;
    const int size_x = nii_rim->nx;
    const int size_y = nii_rim->ny;

    const int nr_voxels = size_z * size_y * size_x;

    const float dX = nii_rim->pixdim[1];
    const float dY = nii_rim->pixdim[2];
    const float dZ = nii_rim->pixdim[3];

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Prepare required nifti images

    nifti_image* fromWM_steps = copy_nifti_header_as_float(nii_rim);
    float* fromWM_steps_data = static_cast<float*>(fromWM_steps->data);
    nifti_image* fromWM_dist = copy_nifti_header_as_float(nii_rim);
    float* fromWM_dist_data = static_cast<float*>(fromWM_dist->data);

    nifti_image* fromGM_steps = copy_nifti_header_as_float(nii_rim);
    float* fromGM_steps_data = static_cast<float*>(fromGM_steps->data);
    nifti_image* fromGM_dist = copy_nifti_header_as_float(nii_rim);
    float* fromGM_dist_data = static_cast<float*>(fromGM_dist->data);

    nifti_image* fromWM_id = copy_nifti_header_as_uint(nii_rim);
    unsigned int* fromWM_id_data = static_cast<unsigned int*>(fromWM_id->data);
    nifti_image* fromGM_id = copy_nifti_header_as_uint(nii_rim);
    unsigned int* fromGM_id_data = static_cast<unsigned int*>(fromGM_id->data);

    nifti_image* layers_equidist  = copy_nifti_header_as_int(nii_rim);
    int* layers_equidist_data = static_cast<int*>(layers_equidist->data);

    nifti_image* nii_columns = copy_nifti_header_as_int(nii_rim);
    int* nii_columns_data = static_cast<int*>(nii_columns->data);

    nifti_image* middle_gm = copy_nifti_header_as_int(nii_rim);
    int* middle_gm_data = static_cast<int*>(middle_gm->data);

    // Setting zero
    for (int i = 0; i != nr_voxels; ++i) {
        *(fromWM_steps_data + i) = 0;
        *(fromWM_id_data + i) = 0;
        *(fromGM_steps_data + i) = 0;
        *(fromGM_id_data + i) = 0;
        *(middle_gm_data + i) = 0;
    }

    // ========================================================================
    // Grow from WM
    // ========================================================================
    cout << "  Start growing from inner GM (WM-facing border)..." << endl;

    // Initialize grow volume
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 2) {  // WM boundary voxels within GM
            *(fromWM_steps_data + i) = 1.;
            *(fromWM_dist_data + i) = 0.;
            *(fromWM_id_data + i) = i;
        } else {
            *(fromWM_steps_data + i) = 0.;
            *(fromWM_dist_data + i) = 0.;
        }
    }

    unsigned int grow_step = 1, voxel_counter = nr_voxels;
    int ix, iy, iz, j;
    float d;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(fromWM_steps_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dX;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dX;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dY;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dY;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dZ;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dZ;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_step + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    save_output_nifti(fin, "fromWM_steps", fromWM_steps, false);
    save_output_nifti(fin, "fromWM_dist", fromWM_dist, false);
    save_output_nifti(fin, "fromWM_id", fromWM_dist, false);

    // ========================================================================
    // Grow from CSF
    // ========================================================================
    cout << "  Start growing from outer GM..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 1) {
            *(fromGM_steps_data + i) = 1.;
            *(fromGM_dist_data + i) = 0.;
            *(fromGM_id_data + i) = i;
        } else {
            *(fromGM_steps_data + i) = 0.;
            *(fromGM_dist_data + i) = 0.;
        }
    }

    grow_step = 1, voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(fromGM_steps_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dX;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dX;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dY;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dY;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dZ;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dZ;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 2) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_step + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    save_output_nifti(fin, "fromGM_steps", fromGM_steps, false);
    save_output_nifti(fin, "fromGM_dist", fromGM_dist, false);
    save_output_nifti(fin, "fromGM_id", fromGM_dist, false);

    // ========================================================================
    // Layers
    // ========================================================================
    cout << "  Doing layers..." << endl;
    float x, y, z, wm_x, wm_y, wm_z, gm_x, gm_y, gm_z, mid_x, mid_y, mid_z;

    // Repurpose data arrays
    for (int i = 0; i != nr_voxels; ++i) {
        *(fromWM_dist_data + i) = 0.;
        *(fromGM_dist_data + i) = 0.;
    }

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(fromWM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(fromGM_id_data + i),
                                               size_x, size_y);

            // Normalize distance
            float dist1, dist2;
            dist1 = dist(x, y, z, wm_x, wm_y, wm_z, dX, dY, dZ);
            dist2 = dist(x, y, z, gm_x, gm_y, gm_z, dX, dY, dZ);
            float norm_dist = dist1 / (dist1 + dist2);

            // Cast distances to integers as number of desired layers
            *(layers_equidist_data + i) = ceil(nr_layers * norm_dist);

            // Middle gray matter
            mid_x = round((gm_x + wm_x) / 2.);
            mid_y = round((gm_y + wm_y) / 2.);
            mid_z = round((gm_z + wm_z) / 2.);
            j = sub2ind_3D(mid_x, mid_y, mid_z, size_x, size_y);
            *(middle_gm_data + j) = 1;

            j = *(fromWM_id_data + i);
            *(fromWM_dist_data + j) += 1;
            j = *(fromGM_id_data + i);
            *(fromGM_dist_data + j) += 1;
        }
    }
    save_output_nifti(fin, "layers_equidist", layers_equidist);
    save_output_nifti(fin, "middle_gm", middle_gm);
    save_output_nifti(fin, "WM_hotspots", fromWM_dist);
    save_output_nifti(fin, "GM_hotspots", fromGM_dist);

    // ========================================================================
    // Columns
    // ========================================================================
    cout << "  Doing columns..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(fromWM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(fromGM_id_data + i),
                                               size_x, size_y);

            // Find middle point of columns
            mid_x = (wm_x + gm_x) / 2.;
            mid_y = (wm_y + gm_y) / 2.;
            mid_z = (wm_z + gm_z) / 2.;

            // Downsample middle point coordinate (makes columns larger)
            mid_x = round(mid_x / column_size) * column_size;
            mid_y = round(mid_y / column_size) * column_size;
            mid_z = round(mid_z / column_size) * column_size;

            j = sub2ind_3D(mid_x, mid_y, mid_z, size_x, size_y);
            *(nii_columns_data + i) = j;
        }
    }
    save_output_nifti(fin, "columns", nii_columns);

    cout << "  Finished." << endl;
    return 0;
}
