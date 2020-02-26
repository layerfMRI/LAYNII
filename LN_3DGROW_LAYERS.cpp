
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
    nifti_image*nim_input = NULL;
    char* fin = NULL;
    int ac, nr_layers = 3;
    if (argc < 2) {
        return show_help();   // typing '-help' is sooo much work
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-nr_layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
            } else {
                nr_layers = atof(argv[ac]);
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
    nim_input = nifti_image_read(fin, 1);
    if (!nim_input) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }
    int16_t *nim_input_data = static_cast<__int16_t*>(nim_input->data);

    log_welcome("LN_3DGROW_LAYERS");
    log_nifti_descriptives(nim_input);

    cout << "  Nr. layers: " << nr_layers << endl;

    // NOTE(Faruk): This is mostly redundant now, probably will take out
    int vinc = 200;

    // Get dimensions of input
    const int size_z = nim_input->nz;
    const int size_x = nim_input->nx;
    const int size_y = nim_input->ny;

    const int nr_voxels = size_z * size_y * size_x;

    const float dX = nim_input->pixdim[1];
    const float dY = nim_input->pixdim[2];
    const float dZ = nim_input->pixdim[3];

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Prepare required nifti images

    nifti_image* fromWM_steps = copy_nifti_header_as_float(nim_input);
    float* fromWM_steps_data = static_cast<float*>(fromWM_steps->data);
    nifti_image* fromWM_dist = copy_nifti_header_as_float(nim_input);
    float* fromWM_dist_data = static_cast<float*>(fromWM_dist->data);

    nifti_image* fromGM_steps = copy_nifti_header_as_float(nim_input);
    float* fromGM_steps_data = static_cast<float*>(fromGM_steps->data);
    nifti_image* fromGM_dist = copy_nifti_header_as_float(nim_input);
    float* fromGM_dist_data = static_cast<float*>(fromGM_dist->data);

    nifti_image* fromWM_id = copy_nifti_header_as_uint(nim_input);
    unsigned int* fromWM_id_data = static_cast<unsigned int*>(fromWM_id->data);
    nifti_image* fromGM_id = copy_nifti_header_as_uint(nim_input);
    unsigned int* fromGM_id_data = static_cast<unsigned int*>(fromGM_id->data);

    nifti_image* layers_equidist  = copy_nifti_header_as_int(nim_input);
    int* layers_equidist_data = static_cast<int*>(layers_equidist->data);

    nifti_image* nii_columns = copy_nifti_header_as_int(nim_input);
    int* nii_columns_data = static_cast<int*>(nii_columns->data);

    // Setting zero
    for (int i = 0; i != nr_voxels; ++i) {
        *(fromWM_steps_data + i) = 0;
        *(fromWM_id_data + i) = 0;
        *(fromGM_steps_data + i) = 0;
        *(fromGM_id_data + i) = 0;
    }

    // ========================================================================
    // Grow from WM
    // ========================================================================
    cout << "  Start growing from inner GM (WM-facing border)..." << endl;

    // Initialize grow volume
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 2) {  // WM boundary voxels within GM
            *(fromWM_steps_data + i) = 1.;
            *(fromWM_dist_data + i) = 1.;
            *(fromWM_id_data + i) = i;
        } else {
            *(fromWM_steps_data + i) = 0.;
            *(fromWM_dist_data + i) = 0.;
        }
    }

    for (int grow_i = 1; grow_i != vinc; grow_i++) {
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(fromWM_steps_data + i) == grow_i) {
                int ix, iy, iz;
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                int j;
                float d;
                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dX;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dX;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dY;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dY;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dZ;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dZ;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xy;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_yz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromWM_dist_data + i) + dia_xyz;
                        if (d < *(fromWM_dist_data + j)
                            || *(fromWM_dist_data + j) == 0) {
                            *(fromWM_dist_data + j) = d;
                            *(fromWM_steps_data + j) = grow_i + 1;
                            *(fromWM_id_data + j) = *(fromWM_id_data + i);
                        }
                    }
                }
            }
        }
    }
    save_output_nifti(fin, "fromWM_steps", fromWM_steps, false);
    save_output_nifti(fin, "fromWM_dist", fromWM_dist, false);
    save_output_nifti(fin, "fromWM_id", fromWM_dist, false);

    // ========================================================================
    // Grow from CSF
    // ========================================================================
    cout << "  Start growing from outer GM..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 1) {
            *(fromGM_steps_data + i) = 1.;
            *(fromGM_dist_data + i) = 1.;
            *(fromGM_id_data + i) = i;
        } else {
            *(fromGM_steps_data + i) = 0.;
            *(fromGM_dist_data + i) = 0.;
        }
    }

    for (int grow_i = 1; grow_i != vinc; grow_i++) {
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(fromGM_steps_data + i) == grow_i) {
                int ix, iy, iz;
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                int j;
                float d;
                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dX;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dX;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dY;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dY;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dZ;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dZ;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xy;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_yz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(fromGM_dist_data + i) + dia_xyz;
                        if (d < *(fromGM_dist_data + j)
                            || *(fromGM_dist_data + j) == 0) {
                            *(fromGM_dist_data + j) = d;
                            *(fromGM_steps_data + j) = grow_i + 1;
                            *(fromGM_id_data + j) = *(fromGM_id_data + i);
                        }
                    }
                }
            }
        }
    }
    save_output_nifti(fin, "fromWM_steps", fromGM_steps, false);
    save_output_nifti(fin, "fromWM_dist", fromGM_dist, false);
    save_output_nifti(fin, "fromWM_id", fromGM_dist, false);

    // ========================================================================
    // Layers
    // ========================================================================
    cout << "  Doing layers..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 3) {
            float x, y, z, wm_x, wm_y, wm_z, gm_x, gm_y, gm_z;
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
        }
    }
    save_output_nifti(fin, "layers_equidist", layers_equidist);

    // ========================================================================
    // Columns
    // ========================================================================
    cout << "  Doing columns..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 3) {
            float x, y, z, wm_x, wm_y, wm_z, gm_x, gm_y, gm_z;
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(fromWM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(fromGM_id_data + i),
                                               size_x, size_y);

            // NOTE(Faruk): Columns, WIP...
            int mid_x = (wm_x + gm_x)/2;
            int mid_y = (wm_y + gm_y)/2;
            int mid_z = (wm_z + gm_z)/2;
            int j = sub2ind_3D(mid_x, mid_y, mid_z, size_x, size_y);

            *(nii_columns_data + i) = round(j);
        }
    }
    save_output_nifti(fin, "columns", nii_columns);

    cout << "  Finished." << endl;
    return 0;
}
