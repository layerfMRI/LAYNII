
// TODO(Faruk): Make columns optional output.
// TODO(Faruk): Put neighbour visits into a function.
// TODO(Faruk): Add debug flag to reduce number of outputs.

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_LAYERS: Cortical gray matter layering.\n"
    "\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "!! BEWARE! WORK IN PROGRESS... USE WITH CAUTION !!\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "\n"
    "Usage:\n"
    "    LN2_LAYERS -rim rim.nii \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help. \n"
    "    -rim        : Specify input dataset.\n"
    "    -nr_layers  : Number of layers. Default is 3.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image*nii1 = NULL;
    char* fin = NULL;
    uint16_t ac, nr_layers = 3;
    float column_size = 1;
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
    nii1 = nifti_image_read(fin, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_3DGROW_LAYERS");
    log_nifti_descriptives(nii1);

    cout << "  Nr. layers: " << nr_layers << endl;

    // Get dimensions of input
    const uint32_t size_z = nii1->nz;
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_rim = copy_nifti_as_int32(nii1);
    int32_t* nii_rim_data = static_cast<int32_t*>(nii_rim->data);

    // Prepare required nifti images
    nifti_image* nii_layers  = copy_nifti_as_int32(nii_rim);
    int32_t* nii_layers_data = static_cast<int32_t*>(nii_layers->data);
    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_layers_data + i) = 0;
    }

    nifti_image* innerGM_step = copy_nifti_as_float32(nii_layers);
    float* innerGM_step_data = static_cast<float*>(innerGM_step->data);
    nifti_image* innerGM_dist = copy_nifti_as_float32(nii_layers);
    float* innerGM_dist_data = static_cast<float*>(innerGM_dist->data);

    nifti_image* outerGM_step = copy_nifti_as_float32(nii_layers);
    float* outerGM_step_data = static_cast<float*>(outerGM_step->data);
    nifti_image* outerGM_dist = copy_nifti_as_float32(nii_layers);
    float* outerGM_dist_data = static_cast<float*>(outerGM_dist->data);

    nifti_image* err_dist = copy_nifti_as_float32(nii_layers);
    float* err_dist_data = static_cast<float*>(err_dist->data);

    nifti_image* innerGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* innerGM_id_data = static_cast<int32_t*>(innerGM_id->data);
    nifti_image* outerGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* outerGM_id_data = static_cast<int32_t*>(outerGM_id->data);

    nifti_image* innerGM_prevstep_id = copy_nifti_as_int32(nii_layers);
    int32_t* innerGM_prevstep_id_data =
        static_cast<int32_t*>(innerGM_prevstep_id->data);
    nifti_image* outerGM_prevstep_id = copy_nifti_as_int32(nii_layers);
    int32_t* outerGM_prevstep_id_data =
        static_cast<int32_t*>(outerGM_prevstep_id->data);
    nifti_image* normdistdiff = copy_nifti_as_float32(nii_layers);
    float* normdistdiff_data = static_cast<float*>(normdistdiff->data);

    nifti_image* nii_columns = copy_nifti_as_int32(nii_layers);
    int32_t* nii_columns_data = static_cast<int32_t*>(nii_columns->data);

    nifti_image* midGM = copy_nifti_as_int32(nii_layers);
    int32_t* midGM_data = static_cast<int32_t*>(midGM->data);
    nifti_image* midGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* midGM_id_data = static_cast<int32_t*>(midGM_id->data);

    nifti_image* hotspots = copy_nifti_as_int32(nii_layers);
    int32_t* hotspots_data = static_cast<int32_t*>(hotspots->data);
    nifti_image* curvature = copy_nifti_as_float32(nii_layers);
    float* curvature_data = static_cast<float*>(curvature->data);
    nifti_image* thickness = copy_nifti_as_float32(nii_layers);
    float* thickness_data = static_cast<float*>(thickness->data);

    // ========================================================================
    // Grow from WM
    // ========================================================================
    cout << "  Start growing from inner GM (WM-facing border)..." << endl;

    // Initialize grow volume
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 2) {  // WM boundary voxels within GM
            *(innerGM_step_data + i) = 1.;
            *(innerGM_dist_data + i) = 0.;
            *(innerGM_id_data + i) = i;
        } else {
            *(innerGM_step_data + i) = 0.;
            *(innerGM_dist_data + i) = 0.;
        }
    }

    uint32_t grow_step = 1, voxel_counter = nr_voxels;
    uint32_t ix, iy, iz, j, k;
    float d;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(innerGM_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dX;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dX;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dY;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dY;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dZ;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dZ;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    save_output_nifti(fin, "innerGM_step", innerGM_step, false);
    save_output_nifti(fin, "innerGM_dist", innerGM_dist, false);
    save_output_nifti(fin, "innerGM_id", innerGM_id, false);

    // ========================================================================
    // Grow from CSF
    // ========================================================================
    cout << "  Start growing from outer GM..." << endl;

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 1) {
            *(outerGM_step_data + i) = 1.;
            *(outerGM_dist_data + i) = 0.;
            *(outerGM_id_data + i) = i;
        } else {
            *(outerGM_step_data + i) = 0.;
            *(outerGM_dist_data + i) = 0.;
        }
    }

    grow_step = 1, voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(outerGM_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dX;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dX;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dY;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dY;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dZ;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dZ;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_z && iz != 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy != size_z && iz != size_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iz != 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix != 0 && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != 0 && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != 0 && iz != 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != 0 && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix != size_x && iy != size_y && iz != size_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    save_output_nifti(fin, "outerGM_step", outerGM_step, false);
    save_output_nifti(fin, "outerGM_dist", outerGM_dist, false);
    save_output_nifti(fin, "outerGM_id", outerGM_id, false);

    // ========================================================================
    // Layers
    // ========================================================================
    cout << "  Start doing layers..." << endl;
    float x, y, z, wm_x, wm_y, wm_z, gm_x, gm_y, gm_z;

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(innerGM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(outerGM_id_data + i),
                                               size_x, size_y);

            // // Normalize distance
            // float dist1 = dist(x, y, z, wm_x, wm_y, wm_z, dX, dY, dZ);
            // float dist2 = dist(x, y, z, gm_x, gm_y, gm_z, dX, dY, dZ);
            // float norm_dist = dist1 / (dist1 + dist2);

            // Normalize distance (completely discrete)
            float dist1 = *(innerGM_dist_data + i);
            float dist2 = *(outerGM_dist_data + i);
            float total_dist = dist1 + dist2;
            *(thickness_data + i) = total_dist;
            float norm_dist = dist1 / total_dist;

            // Difference of normalized distances
            *(normdistdiff_data + i) = (dist1 - dist2) / total_dist;

            // Cast distances to integers as number of desired layers
            if (norm_dist != 0) {
                *(nii_layers_data + i) = ceil(nr_layers * norm_dist);
            } else {
                *(nii_layers_data + i) = 1;
            }

            // NOTE: for debugging purposes
            float dist3 = dist(wm_x, wm_y, wm_z, gm_x, gm_y, gm_z, dX, dY, dZ);
            *(err_dist_data + i) = (dist1 + dist2) - dist3;

            // Count inner and outer GM anchor voxels
            j = *(innerGM_id_data + i);
            *(hotspots_data + j) += 1;
            j = *(outerGM_id_data + i);
            *(hotspots_data + j) -= 1;
        }
    }
    save_output_nifti(fin, "layers_equidist", nii_layers);
    save_output_nifti(fin, "thickness", thickness, false);
    save_output_nifti(fin, "hotspots", hotspots, false);
    save_output_nifti(fin, "disterror", err_dist, false);
    save_output_nifti(fin, "normdistdiff_equidist", normdistdiff, false);

    // ========================================================================
    // Middle gray matter
    // ========================================================================
    cout << "  Start finding middle gray matter..." << endl;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            // Check sign changes in normalized distance differences between
            // neighbouring voxels on a column path (a.k.a. streamline)
            if (*(normdistdiff_data + i) == 0) {
                *(midGM_data + i) = 1;
                *(midGM_id_data + i) = i;
            } else {
                float m = *(normdistdiff_data + i);
                float n;

                // Inner neighbour
                j = *(innerGM_prevstep_id_data + i);
                if (*(nii_rim_data + j) == 3) {
                    n = *(normdistdiff_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (abs(m) < abs(n)) {
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                        } else if (abs(m) > abs(n)) {  // Closer to prev. step
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = j;
                        } else {  // Equal +/- normalized distance
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = i;  // On purpose
                        }
                    }
                }

                // Outer neighbour
                j = *(outerGM_prevstep_id_data + i);
                if (*(nii_rim_data + j) == 3) {
                    n = *(normdistdiff_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (abs(m) < abs(n)) {
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                        } else if (abs(m) > abs(n)) {  // Closer to prev. step
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = j;
                        } else {  // Equal +/- normalized distance
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = i;  // On purpose
                        }
                    }
                }
            }
        }
    }
    save_output_nifti(fin, "midGM", midGM, false);

    // ========================================================================
    // Columns
    // ========================================================================
    cout << "  Start doing columns..." << endl;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            // Approximate curvature measurement per column/streamline
            j = *(innerGM_id_data + i);
            k = *(outerGM_id_data + i);  // These values are negative
            *(curvature_data + i) = *(hotspots_data + j) + *(hotspots_data + k);
            *(curvature_data + i) /=
                *(hotspots_data + j) - *(hotspots_data + k);  // normalize

            // Re-assign mid-GM id based on curvature
            if (*(curvature_data + i) >= 0) {  // Gyrus
                *(nii_columns_data + i) = j;
            } else {
                *(nii_columns_data + i) = k;
            }

            // MiddleGM ids are used to find centroids in the next step
            if (*(midGM_data + i) == 1) {
                *(midGM_id_data + i) = *(nii_columns_data + i);
            }
        }
    }

    // ========================================================================
    // Find Middle Gray Matter centroids
    // ========================================================================
    // NOTE(Faruk): I am a bit sluggish with memory usage here. Might optimize
    // later by switching nifti images to vectors.
    nifti_image* coords_x = copy_nifti_as_float32(nii_rim);
    float* coords_x_data = static_cast<float*>(coords_x->data);
    nifti_image* coords_y = copy_nifti_as_float32(nii_rim);
    float* coords_y_data = static_cast<float*>(coords_y->data);
    nifti_image* coords_z = copy_nifti_as_float32(nii_rim);
    float* coords_z_data = static_cast<float*>(coords_z->data);
    nifti_image* coords_count = copy_nifti_as_int32(nii_rim);
    int32_t* coords_count_data = static_cast<int32_t*>(coords_count->data);
    nifti_image* centroid = copy_nifti_as_int32(nii_rim);
    int32_t* centroid_data = static_cast<int32_t*>(centroid->data);
    nifti_image* midGM_centroid_id = copy_nifti_as_int32(nii_columns);
    int32_t* midGM_centroid_id_data = static_cast<int32_t*>(midGM_centroid_id->data);

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(coords_x_data + i) = 0;
        *(coords_y_data + i) = 0;
        *(coords_z_data + i) = 0;
        *(coords_count_data + i) = 0;
        *(centroid_data + i) = 0;
    }

    // Sum x, y, z coordinates of same-column middle GM voxels
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(midGM_data + i) == 1) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            j = *(midGM_id_data + i);  // used to determine storage voxel
            *(coords_x_data + j) += x;
            *(coords_y_data + j) += y;
            *(coords_z_data + j) += z;
            *(coords_count_data + j) += 1;
        }
    }
    // Divide summed coordinates to find centroid
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(coords_count_data + i) != 0) {
            // Assign centroid id in place of inner/outer border voxel id
            x = floor(*(coords_x_data + i) / *(coords_count_data + i));
            y = floor(*(coords_y_data + i) / *(coords_count_data + i));
            z = floor(*(coords_z_data + i) / *(coords_count_data + i));
            j = sub2ind_3D(x, y, z, size_x, size_y);
            *(centroid_data + i) = j;
        }
    }
    // Map new centroid IDs to columns/streamlines
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(nii_columns_data + i);
            *(midGM_centroid_id_data + i) = *(centroid_data + j);

            if (*(midGM_data + i) == 1) {  // Update Mid GM id
                *(midGM_id_data + i) = *(centroid_data + j);
            }
        }
    }
    save_output_nifti(fin, "midGM_id", midGM_id, false);
    save_output_nifti(fin, "columns", midGM_centroid_id, false);

    // ========================================================================
    // Update curvature along column/streamline based on midGM curvature
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(midGM_centroid_id_data + i);
            *(curvature_data + i) = *(curvature_data + j);
        }
    }
    save_output_nifti(fin, "curvature", curvature, false);

    // ========================================================================
    // Find column centroids
    // ========================================================================
    nifti_image* column_centroid_id = copy_nifti_as_int32(nii_columns);
    int32_t* column_centroid_id_data = static_cast<int32_t*>(column_centroid_id->data);

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(coords_x_data + i) = 0;
        *(coords_y_data + i) = 0;
        *(coords_z_data + i) = 0;
        *(coords_count_data + i) = 0;
        *(centroid_data + i) = 0;
    }

    // Sum x, y, z coordinates of same-column voxels
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            j = *(nii_columns_data + i);  // used to determine storage voxel
            *(coords_x_data + j) += x;
            *(coords_y_data + j) += y;
            *(coords_z_data + j) += z;
            *(coords_count_data + j) += 1;
        }
    }
    // Divide summed coordinates to find centroid
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(coords_count_data + i) != 0) {
            // Assign centroid id in place of inner/outer border voxel id
            x = floor(*(coords_x_data + i) / *(coords_count_data + i));
            y = floor(*(coords_y_data + i) / *(coords_count_data + i));
            z = floor(*(coords_z_data + i) / *(coords_count_data + i));
            j = sub2ind_3D(x, y, z, size_x, size_y);
            *(centroid_data + i) = j;
        }
    }
    // Map new centroid IDs to columns/streamlines
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(nii_columns_data + i);
            *(column_centroid_id_data + i) = *(centroid_data + j);
        }
    }
    // save_output_nifti(fin, "columns2", column_centroid_id, false);

    // ------------------------------------------------------------------------
    // Debug centroids
    nifti_image *debug = copy_nifti_as_int32(nii_rim);
    int32_t *debug_data = static_cast<int32_t*>(debug->data);
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(midGM_centroid_id_data + i);
            *(debug_data + j) = 1;
            j = *(column_centroid_id_data + i);
            *(debug_data + j) = 2;
        } else {
            *(debug_data + i) = 0;
        }
    }
    save_output_nifti(fin, "centroid_debug", debug, false);
    // ------------------------------------------------------------------------

    // ========================================================================
    // Equi-volume
    // ========================================================================
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_rim_data + i) == 3) {
    //         // Normalize distance
    //         float dist1 = *(innerGM_dist_data + i);
    //         float dist2 = *(outerGM_dist_data + i);
    //         float total_dist = (dist1 + dist2);
    //         *(innerGM_dist_data + i) /= total_dist;
    //         *(outerGM_dist_data + i) /= total_dist;
    //     }
    // }
    //
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_rim_data + i) == 3) {
    //         j = *(midGM_centroid_id_data + i);
    //         k = *(column_centroid_id_data + j);
    //
    //         float j1 = *(innerGM_dist_data + j);
    //         float j2 = *(outerGM_dist_data + j);
    //
    //         float k1 = *(innerGM_dist_data + k);
    //         float k2 = *(outerGM_dist_data + k);
    //
    //         float a = j1 / k1;
    //         float b = j2 / k2;
    //
    //         // Normalize distance (completely discrete)
    //         float dist1 = *(innerGM_dist_data + i) * exp(a);
    //         float dist2 = *(outerGM_dist_data + i) * exp(b);
    //         float norm_dist = dist1 / (dist1 + dist2);
    //
    //         // Difference of normalized distances
    //         *(normdistdiff_data + i) = (dist1 - dist2) / (dist1 + dist2);
    //
    //         // Cast distances to integers as number of desired layers
    //         if (norm_dist != 0) {
    //             *(nii_layers_data + i) = ceil(nr_layers * norm_dist);
    //         } else {
    //             *(nii_layers_data + i) = 1;
    //         }
    //     }
    // }
    // save_output_nifti(fin, "layers_equivol", nii_layers);
    // save_output_nifti(fin, "normdistdiff_equivol", normdistdiff, false);

    // ========================================================================
    // Middle gray matter for equi-volume
    // ========================================================================
    // cout << "  Start finding middle gray matter..." << endl;
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     *(midGM_data + i) = 0;
    //     *(midGM_id_data + i) = 0;
    // }
    //
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_rim_data + i) == 3) {
    //         // Check sign changes in normalized distance differences between
    //         // neighbouring voxels on a column path (a.k.a. streamline)
    //         if (*(normdistdiff_data + i) == 0) {
    //             *(midGM_data + i) = 1;
    //             *(midGM_id_data + i) = i;
    //         } else {
    //             float m = *(normdistdiff_data + i);
    //             float n;
    //
    //             // Inner neighbour
    //             j = *(innerGM_prevstep_id_data + i);
    //             if (*(nii_rim_data + j) == 3) {
    //                 n = *(normdistdiff_data + j);
    //                 if (signbit(m) - signbit(n) != 0) {
    //                     if (abs(m) < abs(n)) {
    //                         *(midGM_data + i) = 1;
    //                         *(midGM_id_data + i) = i;
    //                     } else if (abs(m) > abs(n)) {  // Closer to prev. step
    //                         *(midGM_data + j) = 1;
    //                         *(midGM_id_data + j) = j;
    //                     } else {  // Equal +/- normalized distance
    //                         *(midGM_data + i) = 1;
    //                         *(midGM_id_data + i) = i;
    //                         *(midGM_data + j) = 1;
    //                         *(midGM_id_data + j) = i;  // On purpose
    //                     }
    //                 }
    //             }
    //
    //             // Outer neighbour
    //             j = *(outerGM_prevstep_id_data + i);
    //             if (*(nii_rim_data + j) == 3) {
    //                 n = *(normdistdiff_data + j);
    //                 if (signbit(m) - signbit(n) != 0) {
    //                     if (abs(m) < abs(n)) {
    //                         *(midGM_data + i) = 1;
    //                         *(midGM_id_data + i) = i;
    //                     } else if (abs(m) > abs(n)) {  // Closer to prev. step
    //                         *(midGM_data + j) = 1;
    //                         *(midGM_id_data + j) = j;
    //                     } else {  // Equal +/- normalized distance
    //                         *(midGM_data + i) = 1;
    //                         *(midGM_id_data + i) = i;
    //                         *(midGM_data + j) = 1;
    //                         *(midGM_id_data + j) = i;  // On purpose
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // save_output_nifti(fin, "midGM_equivol", midGM, false);

    // ========================================================================
    // TODO(Faruk): Might use bspline weights to smooth curvature maps a bit.
    // TODO(Faruk): Might be better to use step 1 id's to define columns.
    // TODO(Faruk): Mean coordinate of columns might be used in
    // equivolume layering

    cout << "  Finished." << endl;
    return 0;
}
