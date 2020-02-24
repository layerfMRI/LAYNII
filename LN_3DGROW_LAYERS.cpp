
// TODO(@Renzo): Check the vincinity loops. the "+1" are missing and
// the z_dimensions is limited by REAN dimension????


#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_3DGROW_LAYERS: Short exmample of cortical gray matter layering.\n"
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
    char* fin = NULL, *fout = NULL;
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
    float dZ = nim_input->pixdim[3];

    // Get access to data of nim_input
    if (nim_input->datatype != 4) {
        // nim_input->datatype = NIFTI_TYPE_INT16;
        cout << "  !!!WRONG DATATYPE!!!" << endl;
    }

    // ========================================================================
    // Fix datatype issues

    nifti_image* growfromWM0 = copy_nifti_header_as_float(nim_input);
    float* growfromWM0_data = static_cast<float*>(growfromWM0->data);

    nifti_image* growfromWM1 = copy_nifti_header_as_float(nim_input);
    float* growfromWM1_data = static_cast<float*>(growfromWM1->data);

    nifti_image* growfromGM0 = copy_nifti_header_as_float(nim_input);
    float* growfromGM0_data = static_cast<float*>(growfromGM0->data);

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
    float x2g = 0., y2g = 0., z2g = 0.;
    float x3g = 0., y3g = 0., z3g = 0.;

    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float angle(float a, float b, float c);

    // Reduce mask to contain only Areas close to the curface.
    cout << "  Select GM regions..." << endl;

    // This is the distance from every voxel that the algorithm is applied on.
    // Just to make it faster and not loop over all voxels.
    int vinc = 80;

    float dist_i = 0.;
    float dist_min = 0.;
    float dist_min1 = 0.;
    float dist_min2 = 0.;
    float dist_min3 = 0.;
    float dist_max = 0.;
    float dist_p1 = 0.;

    int nr_layers = 20;

    cout << "  Start growing from WM..." << endl;

    // Setting zero
    for (int i = 0; i != nr_voxels; ++i) {
        *(growfromWM0_data + i) = 0;
        *(growfromWM1_data + i) = 0;
        *(growfromGM0_data + i) = 0;
        *(growfromGM1_data + i) = 0;
    }

    /////////////////////////
    // Closing 3D surfaces //
    /////////////////////////

    // cout << "  Closing surfaces..." << endl;
    //
    // int isWMb, isCSFb, isb;
    //
    // for (int iterations = 0; iterations< 100; ++iterations) {
    //     for (int iz=0; iz<size_z; ++iz) {
    //         for (int iy=0; iy<size_x; ++iy) {
    //             for (int ix=0; ix<size_y-0; ++ix) {
    //                 if (*(nim_input_data + VOXEL_ID_3D) == 0) {
    //                     isWMb=0;
    //                     isCSFb=0;
    //                     isb=0;
    //                     for (int iz_i=max(0,iz-2); iz_i<min(size_z,iz+2); ++iz_i) {
    //                         for (int iy_i=max(0,iy-2); iy_i<minxy * izn(size_x,iy+2); ++iy_i) {
    //                             for (int ix_i=max(0,ix-2); ix_i<min(size_y,ix+2); ++ix_i) {
    //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 3) isb = 1;
    //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 1) isCSFb = 1;
    //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 2) isWMb = 1;
    //                             }
    //                         }
    //                     }
    //                     if (isb == 1 && isCSFb   == 1) *(nim_input_data + VOXEL_ID_3D) = 1;
    //                     if (isb == 1 &&  isWMb   == 1) *(nim_input_data + VOXEL_ID_3D) = 2;
    //                     if (isCSFb == 1 && isWMb == 1) {
    //                         *(nim_input_data + VOXEL_ID_3D) = 0;
    //                         cout << " THIS BRAIN IS WEIRD." << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    //////////////////
    // Grow from WM //
    //////////////////
    int grow_vinc = 2;

    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nim_input_data + i) == 2) {
            int x, y, z;
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y, size_z);
            *(growfromWM0_data + i) = 1.;
            *(WMkoordx1_data + i) = x;
            *(WMkoordy1_data + i) = y;
            *(WMkoordz1_data + i) = z;
        }
    }

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (*(nim_input_data + VOXEL_ID_3D) == 3 && *(growfromWM0_data + VOXEL_ID_3D) == 0) {
                        // cout << "  True   " << *(growfromWM0_data + nxy*islice + nx * ix + iy)<< endl;
                        for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc, size_z); ++iz_i) {
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc, size_y); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc, size_x); ++ix_i) {
                                    if (*(growfromWM0_data + nxy * iz_i + nx * iy_i + ix_i) == (float)grow_i) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        if (dist_i < dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            z1g = iz_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                        }
                        if (dist_min2 < 1.4) {
                            // distDebug(0,islice,iy,ix) = dist_min2;
                            *(growfromWM0_data + VOXEL_ID_3D) = (float)grow_i + 1;
                            *(WMkoordx1_data + VOXEL_ID_3D) = *(WMkoordx1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                            *(WMkoordy1_data + VOXEL_ID_3D) = *(WMkoordy1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                            *(WMkoordz1_data + VOXEL_ID_3D) = *(WMkoordz1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                        }
                    }
                }
            }
        }
    }

    ///////////////////
    // Grow from CSF //
    ///////////////////
    cout << "  Start growing from CSF ..." << endl;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                if (*(nim_input_data + VOXEL_ID_3D) == 1) {
                    *(growfromGM0_data + VOXEL_ID_3D) = 1.;
                    *(GMkoordx1_data + VOXEL_ID_3D) = ix;
                    *(GMkoordy1_data + VOXEL_ID_3D) = iy;
                    *(GMkoordz1_data + VOXEL_ID_3D) = iz;
                }
            }
        }
    }

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0.;
                    y1g = 0;
                    z1g = 0;
                    if (*(nim_input_data + VOXEL_ID_3D) == 3 && *(growfromGM0_data + VOXEL_ID_3D) == 0) {
                        for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc, size_z); ++iz_i) {
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc, size_y); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc, size_x); ++ix_i) {
                                    if (*(growfromGM0_data + nxy * iz_i + nx * iy_i + ix_i) == (float)grow_i) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        if (dist_i < dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            z1g = iz_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                        }
                        if (dist_min2 < 1.4) {
                            *(growfromGM0_data + VOXEL_ID_3D) = (float)grow_i + 1;
                            *(GMkoordx1_data + VOXEL_ID_3D) = *(GMkoordx1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                            *(GMkoordy1_data + VOXEL_ID_3D) = *(GMkoordy1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                            *(GMkoordz1_data + VOXEL_ID_3D) = *(GMkoordz1_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                        }
                    }
                }
            }
        }
    }

    ////////////////////////////////////////////////////
    // Wabble across neigbouring voxles of closest WM //
    // to account for Pythagoras errors               //
    ////////////////////////////////////////////////////
    cout << "  Correct for Pythagoras error..." << endl;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                *(growfromWM1_data + VOXEL_ID_3D) = *(growfromWM0_data + VOXEL_ID_3D);
                *(WMkoordx2_data + VOXEL_ID_3D) = *(WMkoordx1_data + VOXEL_ID_3D);
                *(WMkoordy2_data + VOXEL_ID_3D) = *(WMkoordy1_data + VOXEL_ID_3D);
                *(WMkoordz2_data + VOXEL_ID_3D) = *(WMkoordz1_data + VOXEL_ID_3D);
            }
        }
    }

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    if (*(WMkoordy1_data + VOXEL_ID_3D) != 0) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        z1g = 0;

                        for (int iz_i = max(0, (*(WMkoordz2_data + VOXEL_ID_3D)) - grow_vinc); iz_i < min((*(WMkoordz2_data + VOXEL_ID_3D)) + grow_vinc, size_z); ++iz_i) {
                            for (int iy_i = max(0, (*(WMkoordy2_data + VOXEL_ID_3D)) - grow_vinc); iy_i < min((*(WMkoordy2_data + VOXEL_ID_3D)) + grow_vinc, size_y); ++iy_i) {
                                for (int ix_i = max(0, (*(WMkoordx2_data + VOXEL_ID_3D)) - grow_vinc); ix_i < min((*(WMkoordx2_data + VOXEL_ID_3D)) + grow_vinc, size_x); ++ix_i) {
                                    if (*(nim_input_data + nxy * iz_i + nx * iy_i + ix_i) == 2) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        if (dist_i < dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            z1g = iz_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                        }
                        *(growfromWM1_data + VOXEL_ID_3D) = dist((float)ix, (float)iy, (float)iz, (float)x1g, (float)y1g, (float)z1g, dX, dY, dZ);
                        *(WMkoordx2_data + VOXEL_ID_3D) = *(WMkoordx2_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                        *(WMkoordy2_data + VOXEL_ID_3D) = *(WMkoordy2_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                        *(WMkoordz2_data + VOXEL_ID_3D) = *(WMkoordz2_data + nxy * (int)z1g + nx * (int)y1g + (int)x1g);
                    }
                }
            }
        }
    }
    cout << "  Running until stage 1..." << endl;

    //////////////////////////////////////////////////////
    // Wabble accross neigbouring voexles of closest GM //
    // to account for Pythagoras errors                 //
    //////////////////////////////////////////////////////
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                *(growfromGM1_data + VOXEL_ID_3D) = *(growfromGM0_data + VOXEL_ID_3D);
                *(GMkoordx2_data + VOXEL_ID_3D) = *(GMkoordx1_data + VOXEL_ID_3D);
                *(GMkoordy2_data + VOXEL_ID_3D) = *(GMkoordy1_data + VOXEL_ID_3D);
                *(GMkoordz2_data + VOXEL_ID_3D) = *(GMkoordz1_data + VOXEL_ID_3D);
            }
        }
    }
    cout << "  Running until stage 2..." << endl;

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    if (*(GMkoordy1_data + VOXEL_ID_3D) != 0) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        z1g = 0;
                        for (int iz_i = max(0, (*(GMkoordz2_data + VOXEL_ID_3D)) - grow_vinc); iz_i < min((*(GMkoordz2_data + VOXEL_ID_3D)) + grow_vinc, size_z); ++iz_i) {
                            for (int iy_i = max(0, (*(GMkoordy2_data + VOXEL_ID_3D)) - grow_vinc); iy_i < min((*(GMkoordy2_data + VOXEL_ID_3D)) + grow_vinc, size_y); ++iy_i) {
                                for (int ix_i = max(0, (*(GMkoordx2_data + VOXEL_ID_3D)) - grow_vinc); ix_i < min((*(GMkoordx2_data + VOXEL_ID_3D)) + grow_vinc, size_x); ++ix_i) {
                                    if (*(nim_input_data + nxy * iz_i + nx * iy_i + ix_i)  == 1) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        if (dist_i < dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            z1g = iz_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                        }
                        *(growfromGM1_data + VOXEL_ID_3D) = dist((float)ix, (float)iy, (float)iz, (float)x1g, (float)y1g, (float)z1g, dX, dY, dZ);
                        *(GMkoordx2_data + VOXEL_ID_3D) = *(GMkoordx2_data + nxy*(int)z1g + nx * (int)y1g + (int)x1g);
                        *(GMkoordy2_data + VOXEL_ID_3D) = *(GMkoordy2_data + nxy*(int)z1g + nx * (int)y1g + (int)x1g);
                        *(GMkoordz2_data + VOXEL_ID_3D) = *(GMkoordz2_data + nxy*(int)z1g + nx * (int)y1g + (int)x1g);
                    }
                }
            }
        }
    }
    cout << "  Running until stage 3..." << endl;

    int GMK2_i, GMKz2_i, GMK3_i, WMK2_i, WMKz2_i, WMK3_i;
    float GMK2_f, GMKz2_f, GMK3_f, WMK2_f, WMKz2_f, WMK3_f, ix_f, iy_f, iz_f;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                if (*(nim_input_data + VOXEL_ID_3D) == 3) {
                    // equi_dist_layers(0,iz,iy,ix) = 19 * (1 - dist((float)ix, (float)iy, (float)GMkoord(2, iz, iy, ix), (float)GMkoord(3, iz, iy, ix)) / (dist((float)ix,(float)iy,(float)GMkoord(2, iz, iy, ix),(float)GMkoord(3, iz, iy, ix)) + dist((float)ix,(float)iy,(float)WMkoord(2,iz,iy,ix), (float)WMkoord(3, iz, iy, ix)))) + 2;
                    GMK2_i = *(GMkoordx2_data + VOXEL_ID_3D);
                    GMK3_i = *(GMkoordy2_data + VOXEL_ID_3D);
                    GMKz2_i = *(GMkoordz2_data+ VOXEL_ID_3D);

                    WMK2_i = *(WMkoordx2_data + VOXEL_ID_3D);
                    WMK3_i = *(WMkoordy2_data + VOXEL_ID_3D);
                    WMKz2_i = *(WMkoordz2_data + VOXEL_ID_3D);

                    GMK2_f = (float)GMK2_i;
                    GMK3_f = (float)GMK3_i;
                    GMKz2_f = (float)GMKz2_i;

                    WMK2_f = (float)WMK2_i;
                    WMK3_f = (float)WMK3_i;
                    WMKz2_f = (float)WMKz2_i;

                    ix_f = (float)ix;
                    iy_f = (float)iy;
                    iz_f = (float)iz;

                    // cout << " rix_f,iy_f,GMK2_f,GMK3_f " << "   "  << ix_f << "   "  <<iy_f << "   "  <<GMK2_f << "   "  <<*(GMkoordx2_data + VOXEL_ID_3D)<< endl;
                    *(equi_dist_layers_data + VOXEL_ID_3D) = 19 * (1 - dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ) / (dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ) + dist((float)ix, (float)iy, (float)iz, WMK2_f, WMK3_f, WMKz2_f, dX, dY, dZ))) + 2;
                    //*(equi_dist_layers_data + VOXEL_ID_3D) =  100 * dist(ix_f,iy_f,GMK2_f,GMK3_f);
                }
            }
        }
    }
    cout << "  Running until stage 4..." << endl;

    // Cleaning negative layers and layers ov more than 20
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                if (*(nim_input_data + VOXEL_ID_3D) == 1
                    && *(equi_dist_layers_data + VOXEL_ID_3D) == 0) {
                    *(equi_dist_layers_data + VOXEL_ID_3D) = 21;
                }
                if (*(nim_input_data + VOXEL_ID_3D) == 2
                    && *(equi_dist_layers_data + VOXEL_ID_3D) == 0) {
                    *(equi_dist_layers_data + VOXEL_ID_3D) = 1;
                }
            }
        }
    }
    cout << "  Running until stage 4.5..." << endl;

    // equi_dist_layers.autowrite("equi_dist_layers.nii", wopts, &prot);
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
