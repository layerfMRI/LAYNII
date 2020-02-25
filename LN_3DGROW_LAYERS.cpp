
#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_3DGROW_LAYERS: Cortical gray matter layering.\n"
    "\n"
    "Usage:\n"
    "    LN_3DGROW_LAYERS -rim rim.nii \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help. \n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -rim  border        : Specify input dataset.\n"
    "\n"
    "Notes:\n"
    "    - rim.nii file always needs to be in datatype INT16. \n"
    "    - This is 3D. Hence rim.nii file should be dmsmooth in all \n"
    "      three dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image*nim_input = NULL;
    char* fin = NULL;
    int ac, disp_float_eg = 0;
    if (argc < 2) {
        return show_help();   // typing '-help' is sooo much work
    }
    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-disp_float_example")) {
            disp_float_eg = 1;
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // no string copy, just pointer assignment
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

    // Get dimensions of input
    const int size_z = nim_input->nz;
    const int size_x = nim_input->nx;
    const int size_y = nim_input->ny;
    // const int nrep = nim_input->nt;
    const int nx = nim_input->nx;
    const int nxy = nim_input->nx * nim_input->ny;
    // const int nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    const int nr_voxels = size_z * size_y * size_x;

    const float dX = nim_input->pixdim[1];
    const float dY = nim_input->pixdim[2];
    const float dZ = nim_input->pixdim[3];

    const float min_dim = min(min(dX, dY), dZ);

    // Get access to data of nim_input
    if (nim_input->datatype != 4) {
        // nim_input->datatype = NIFTI_TYPE_INT16;
        cout << "  !!!WRONG DATATYPE!!!" << endl;
    }

    // ========================================================================
    // Fix datatype issues

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

    nifti_image* WMkoordx1 = copy_nifti_header_as_int(nim_input);
    nifti_image* WMkoordy1 = copy_nifti_header_as_int(nim_input);
    nifti_image* WMkoordz1 = copy_nifti_header_as_int(nim_input);
    nifti_image* WMkoordx2 = copy_nifti_header_as_int(nim_input);
    nifti_image* WMkoordy2 = copy_nifti_header_as_int(nim_input);
    nifti_image* WMkoordz2 = copy_nifti_header_as_int(nim_input);

    int* WMkoordx1_data = static_cast<int*>(WMkoordx1->data);
    int* WMkoordy1_data = static_cast<int*>(WMkoordy1->data);
    int* WMkoordz1_data = static_cast<int*>(WMkoordz1->data);
    int* WMkoordx2_data = static_cast<int*>(WMkoordx2->data);
    int* WMkoordy2_data = static_cast<int*>(WMkoordy2->data);
    int* WMkoordz2_data = static_cast<int*>(WMkoordz2->data);

    nifti_image* GMkoordx1 = copy_nifti_header_as_int(nim_input);
    nifti_image* GMkoordy1 = copy_nifti_header_as_int(nim_input);
    nifti_image* GMkoordz1 = copy_nifti_header_as_int(nim_input);
    nifti_image* GMkoordx2 = copy_nifti_header_as_int(nim_input);
    nifti_image* GMkoordy2 = copy_nifti_header_as_int(nim_input);
    nifti_image* GMkoordz2 = copy_nifti_header_as_int(nim_input);

    int* GMkoordx1_data = static_cast<int*>(GMkoordx1->data);
    int* GMkoordy1_data = static_cast<int*>(GMkoordy1->data);
    int* GMkoordz1_data = static_cast<int*>(GMkoordz1->data);
    int* GMkoordx2_data = static_cast<int*>(GMkoordx2->data);
    int* GMkoordy2_data = static_cast<int*>(GMkoordy2->data);
    int* GMkoordz2_data = static_cast<int*>(GMkoordz2->data);

    // ========================================================================

    nifti_image* equi_dist_layers  = copy_nifti_header_as_int(nim_input);
    int* equi_dist_layers_data = static_cast<int*>(equi_dist_layers->data);

    // Coordinates
    float x1g = 0., y1g = 0., z1g = 0.;

    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float angle(float a, float b, float c);

    // Reduce mask to contain only Areas close to the curface.
    cout << "  Select GM regions..." << endl;

    // This is the distance from every voxel that the algorithm is applied on.
    // Just to make it faster and not loop over all voxels.
    int vinc = 200;

    float dist_i = 0.;
    float dist_min2 = 0.;
    float dist_min3 = 0.;
    float dist_max = 0.;

    int nr_layers = 20;
    int grow_vinc = 2;


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
                // 1-jump neighbours
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
            }
        }
    }
    if (nifti_set_filenames(fromWM_steps, "fromWM_steps.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromWM_steps);

    if (nifti_set_filenames(fromWM_dist, "fromWM_dist.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromWM_dist);

    if (nifti_set_filenames(fromWM_id, "fromWM_id.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromWM_id);

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
                // 1-jump neighbours
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
            }
        }
    }
    if (nifti_set_filenames(fromGM_steps, "fromGM_steps.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromGM_steps);

    if (nifti_set_filenames(fromGM_dist, "fromGM_dist.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromGM_dist);

    if (nifti_set_filenames(fromGM_id, "fromGM_id.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(fromGM_id);
    // ========================================================================
    // ========================================================================
    // ========================================================================



    // ========================================================================
    // ========================================================================
    // ========================================================================
    cout << "  Running until stage 3..." << endl;

    int GMK2_i, GMKz2_i, GMK3_i, WMK2_i, WMKz2_i, WMK3_i;
    float GMK2_f, GMKz2_f, GMK3_f, WMK2_f, WMKz2_f, WMK3_f;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 3) {
            GMK2_i  = *(GMkoordx2_data + i);
            GMK3_i  = *(GMkoordy2_data + i);
            GMKz2_i = *(GMkoordz2_data+ i);

            WMK2_i  = *(WMkoordx2_data + i);
            WMK3_i  = *(WMkoordy2_data + i);
            WMKz2_i = *(WMkoordz2_data + i);

            GMK2_f  = static_cast<float>(GMK2_i);
            GMK3_f  = static_cast<float>(GMK3_i);
            GMKz2_f = static_cast<float>(GMKz2_i);

            WMK2_f  = static_cast<float>(WMK2_i);
            WMK3_f  = static_cast<float>(WMK3_i);
            WMKz2_f = static_cast<float>(WMKz2_i);

            int ix, iy, iz;
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);


            float dist1, dist2, dist3;
            dist1 = dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ);
            dist2 = dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ);
            dist3 = dist((float)ix, (float)iy, (float)iz, WMK2_f, WMK3_f, WMKz2_f, dX, dY, dZ);

            *(equi_dist_layers_data + i) = 19 * (1 - dist1 / (dist2 + dist3)) + 2;
        }
    }
    cout << "  Running until stage 4..." << endl;

    // Cleaning negative layers and layers of more than 20
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 1 && *(equi_dist_layers_data + i) == 0) {
            *(equi_dist_layers_data + i) = 21;
        }
        if (*(nim_input_data + i) == 2 && *(equi_dist_layers_data + i) == 0) {
            *(equi_dist_layers_data + i) = 1;
        }
    }
    cout << "  Running until stage 5..." << endl;

    // Output file name
    const char* fout_4 = "equi_dist_layers.nii";
    if (nifti_set_filenames(equi_dist_layers, fout_4, 1, 1)) {
        return 1;
    }
    nifti_image_write(equi_dist_layers);

    const char* fout_5 = "rim_closed.nii";
    if (nifti_set_filenames(nim_input, fout_5, 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_input);

    // const char* fout_6 = "kootrGMx.nii";
    // if (nifti_set_filenames(GMkoordx1, fout_6 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(GMkoordx1);

    // const char* fout_7 = "kootrGMz.nii";
    // if (nifti_set_filenames(GMkoordz1, fout_7 , 1, 1)) {
    //     return 1;
    //     }
    // nifti_image_write(GMkoordz1);

    // koord.autowrite("koordinaten.nii", wopts, &prot);
    cout << "  Finished." << endl;
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt(pow((x1 - x2) * dX, 2) + pow((y1 - y2) * dY, 2)
                + pow((z1 - z2) * dZ, 2));
}

float angle(float a, float b, float c) {
    if (a * a + b * b - c * c <= 0) {
        return 3.141592;
    } else {
        return acos((a * a + b * b - c * c) / (2. * a * b));
    }
}
