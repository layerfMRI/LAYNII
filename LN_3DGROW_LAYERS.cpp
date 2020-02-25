
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

    nifti_image* stepsfromWM0 = copy_nifti_header_as_float(nim_input);
    float* stepsfromWM0_data = static_cast<float*>(stepsfromWM0->data);
    nifti_image* distfromWM0 = copy_nifti_header_as_float(nim_input);
    float* distfromWM0_data = static_cast<float*>(distfromWM0->data);

    nifti_image* growfromWM1 = copy_nifti_header_as_float(nim_input);
    float* growfromWM1_data = static_cast<float*>(growfromWM1->data);

    nifti_image* stepsfromGM0 = copy_nifti_header_as_float(nim_input);
    float* stepsfromGM0_data = static_cast<float*>(stepsfromGM0->data);
    nifti_image* distfromGM0 = copy_nifti_header_as_float(nim_input);
    float* distfromGM0_data = static_cast<float*>(distfromGM0->data);

    nifti_image* growfromGM1 = copy_nifti_header_as_float(nim_input);
    float* growfromGM1_data = static_cast<float*>(growfromGM1->data);

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

    cout << "  Start growing from WM..." << endl;

    // Setting zero
    for (int i = 0; i != nr_voxels; ++i) {
        *(stepsfromWM0_data + i) = 0;
        *(growfromWM1_data + i) = 0;
        *(stepsfromGM0_data + i) = 0;
        *(growfromGM1_data + i) = 0;
    }

    // ========================================================================
    // Grow from WM
    // ========================================================================
    int grow_vinc = 2;

    // Initialize grow volume
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 2) {  // WM boundary voxels within GM
            *(stepsfromWM0_data + i) = 1.;
            *(distfromWM0_data + i) = 1.;
        } else {
            *(stepsfromWM0_data + i) = 0.;
            *(distfromWM0_data + i) = 0.;
        }
    }

    for (int grow_i = 1; grow_i != vinc; grow_i++) {
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(stepsfromWM0_data + i) == grow_i) {
                int ix, iy, iz;
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                int j;
                float d;
                // 1-jump neighbours
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dX;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dX;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dY;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dY;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dZ;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromWM0_data + i) + dZ;
                        if (d < *(distfromWM0_data + j)
                            || *(distfromWM0_data + j) == 0) {
                            *(distfromWM0_data + j) = d;
                            *(stepsfromWM0_data + j) = grow_i + 1;
                        }
                    }
                }
            }
        }
    }
    if (nifti_set_filenames(stepsfromWM0, "stepsfromWM0.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(stepsfromWM0);

    if (nifti_set_filenames(distfromWM0, "distfromWM0.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(distfromWM0);

    // ========================================================================
    // Grow from CSF
    // ========================================================================
    cout << "  Start growing from CSF ..." << endl;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 1) {
            *(stepsfromGM0_data + i) = 1.;
            *(distfromGM0_data + i) = 1.;
        } else {
            *(stepsfromGM0_data + i) = 0.;
            *(distfromGM0_data + i) = 0.;
        }
    }

    for (int grow_i = 1; grow_i != vinc; grow_i++) {
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(stepsfromGM0_data + i) == grow_i) {
                int ix, iy, iz;
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                int j;
                float d;
                // 1-jump neighbours
                if (ix != 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dX;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (ix != size_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dX;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iy != 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dY;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iy != size_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dY;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iz != 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dZ;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
                if (iz != size_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nim_input_data + j) == 3) {
                        d = *(distfromGM0_data + i) + dZ;
                        if (d < *(distfromGM0_data + j)
                            || *(distfromGM0_data + j) == 0) {
                            *(distfromGM0_data + j) = d;
                            *(stepsfromGM0_data + j) = grow_i + 1;
                        }
                    }
                }
            }
        }
    }
    if (nifti_set_filenames(stepsfromGM0, "stepsfromGM0.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(stepsfromGM0);

    if (nifti_set_filenames(distfromGM0, "distfromGM0.nii", 1, 1)) {
        return 1;
    }
    nifti_image_write(distfromGM0);
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
            GMK2_i = *(GMkoordx2_data + i);
            GMK3_i = *(GMkoordy2_data + i);
            GMKz2_i = *(GMkoordz2_data+ i);

            WMK2_i = *(WMkoordx2_data + i);
            WMK3_i = *(WMkoordy2_data + i);
            WMKz2_i = *(WMkoordz2_data + i);

            GMK2_f = static_cast<float>(GMK2_i);
            GMK3_f = static_cast<float>(GMK3_i);
            GMKz2_f = static_cast<float>(GMKz2_i);

            WMK2_f = static_cast<float>(WMK2_i);
            WMK3_f = static_cast<float>(WMK3_i);
            WMKz2_f = static_cast<float>(WMKz2_i);

            int ix, iy, iz;
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

            *(equi_dist_layers_data + i) = 19 * (1 - dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ) / (dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ) + dist((float)ix, (float)iy, (float)iz, WMK2_f, WMK3_f, WMKz2_f, dX, dY, dZ))) + 2;
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
