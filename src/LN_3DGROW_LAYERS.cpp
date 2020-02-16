
// TODO(@Renzo): Check the vincinity loops. the "+1" are missing and
// the z_dimensions is limited by REAN dimension????

#include "./common.h"
#include "./utils.h"

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

    // Get dimensions of input
    const int sizeSlice = nim_input->nz;
    const int sizePhase = nim_input->nx;
    const int sizeRead = nim_input->ny;
    // const int nrep = nim_input->nt;
    const int nx = nim_input->nx;
    const int nxy = nim_input->nx * nim_input->ny;
    // const int nxyz = nim_input->nx * nim_input->ny * nim_input->nz;
    const float dX = nim_input->pixdim[1];
    const float dY = nim_input->pixdim[2];
    const float dZ = nim_input->pixdim[3];

    log_nifti_descriptives(nim_input);

    // if (!fout) { fprintf(stderr, "-- no output requested \n"); return 0; }
    // assign nifti_image fname/iname pair, based on output filename
    //   (request to 'check' image and 'set_byte_order' here)
    // if (nifti_set_filenames(nim_input, fout, 1, 1)) return 1;

    // Get access to data of nim_input
    if (nim_input->datatype != 4) {
        // nim_input->datatype = NIFTI_TYPE_INT16;
        cout << "  !!!WRONG DATATYPE!!!" << endl;
    }

    cout << "  Datatype 1 " << nim_input->datatype << endl;

    // nim_input->datatype = 8;
    // nim_input->nbyper = sizeof(int);
    // nim_input->data = calloc(nim_input->nvox, nim_input->nbyper);

    // float* nim_input_data = (float*) nim_input->data;
    // int* nim_input_data = (int*) nim_input->data;
    // short *nim_input_data = (short *) nim_input->data;
    // int16_t *nim_input_data = (__int16_t *) nim_input->data;
    int16_t *nim_input_data = (__int16_t *) nim_input->data;

    cout << "  Datatype 2 " << nim_input->datatype << endl;

    // nifti_brick_list   NB_orig, NB_select;
    // nifti_image      * nim_orig, * nim_select;
    // const int64_t                blist[5] = { 0, 0, 0, 0, 0 };

    // nim_orig   = nifti_image_read_bricks("rim.nii", 0, NULL,  &NB_orig);
    // nim_select = nifti_image_read_bricks("rim.nii", 5, blist, &NB_select);
    // update_nifti_image_for_brick_list(nim_orig,  1);

    // allocating an additional nii
    // const int64_t blist_2[1] = { 0};
    // nifti_brick_list NB_select_2;
    // nifti_image* growfromWM = nifti_image_read_bricks("rim.nii", 0, NULL, &NB_select_2);

    nifti_image* growfromWM0 = nifti_copy_nim_info(nim_input);
    growfromWM0->datatype = NIFTI_TYPE_FLOAT32;
    growfromWM0->nbyper = sizeof(float);
    growfromWM0->data = calloc(growfromWM0->nvox, growfromWM0->nbyper);
    float* growfromWM0_data = (float*) growfromWM0->data;

    nifti_image* growfromWM1 = nifti_copy_nim_info(nim_input);
    growfromWM1->datatype = NIFTI_TYPE_FLOAT32;
    growfromWM1->nbyper = sizeof(float);
    growfromWM1->data = calloc(growfromWM1->nvox, growfromWM1->nbyper);
    float* growfromWM1_data = (float*) growfromWM1->data;

    // nifti_image* growfromWM1 = nifti_image_read(fin, 1);
    // float* growfromWM1_data = (float*) growfromWM1->data;
    // nifti_image* growfromWM0 = nifti_image_read(fin, 1);
    // float* growfromWM0_data = (float*) growfromWM0->data;
    // growfromWM->dim[4] = 1;
    // nifti_update_dims_from_array(growfromWM);  // changing according sizes nt etc.

    nifti_image* WMkoordx1 = nifti_copy_nim_info(nim_input);
    nifti_image* WMkoordy1 = nifti_copy_nim_info(nim_input);
    nifti_image* WMkoordz1 = nifti_copy_nim_info(nim_input);
    nifti_image* WMkoordx2 = nifti_copy_nim_info(nim_input);
    nifti_image* WMkoordy2 = nifti_copy_nim_info(nim_input);
    nifti_image* WMkoordz2 = nifti_copy_nim_info(nim_input);

    WMkoordx1->datatype = NIFTI_TYPE_INT32;
    WMkoordy1->datatype = NIFTI_TYPE_INT32;
    WMkoordz1->datatype = NIFTI_TYPE_INT32;
    WMkoordx2->datatype = NIFTI_TYPE_INT32;
    WMkoordy2->datatype = NIFTI_TYPE_INT32;
    WMkoordz2->datatype = NIFTI_TYPE_INT32;

    WMkoordx1->nbyper = sizeof(int);
    WMkoordy1->nbyper = sizeof(int);
    WMkoordz1->nbyper = sizeof(int);
    WMkoordx2->nbyper = sizeof(int);
    WMkoordy2->nbyper = sizeof(int);
    WMkoordz2->nbyper = sizeof(int);

    WMkoordx1->data = calloc(WMkoordx1->nvox, WMkoordx1->nbyper);
    WMkoordy1->data = calloc(WMkoordy1->nvox, WMkoordy1->nbyper);
    WMkoordz1->data = calloc(WMkoordz1->nvox, WMkoordz1->nbyper);
    WMkoordx2->data = calloc(WMkoordx2->nvox, WMkoordx2->nbyper);
    WMkoordy2->data = calloc(WMkoordy2->nvox, WMkoordy2->nbyper);
    WMkoordz2->data = calloc(WMkoordz2->nvox, WMkoordz2->nbyper);

    int* WMkoordx1_data = (int*) WMkoordx1->data;
    int* WMkoordy1_data = (int*) WMkoordy1->data;
    int* WMkoordz1_data = (int*) WMkoordz1->data;
    int* WMkoordx2_data = (int*) WMkoordx2->data;
    int* WMkoordy2_data = (int*) WMkoordy2->data;
    int* WMkoordz2_data = (int*) WMkoordz2->data;

    nifti_image* growfromGM0  = nifti_copy_nim_info(nim_input);
    growfromGM0->datatype = NIFTI_TYPE_FLOAT32;
    growfromGM0->nbyper = sizeof(float);
    growfromGM0->data = calloc(growfromGM0->nvox, growfromGM0->nbyper);
    float* growfromGM0_data = (float*) growfromGM0->data;

    nifti_image* growfromGM1  = nifti_copy_nim_info(nim_input);
    growfromGM1->datatype = NIFTI_TYPE_FLOAT32;
    growfromGM1->nbyper = sizeof(float);
    growfromGM1->data = calloc(growfromGM1->nvox, growfromGM1->nbyper);
    float* growfromGM1_data = (float*) growfromGM1->data;

    // nifti_image* growfromGM1 = nifti_image_read(fin, 1);
    // float* growfromGM1_data = (float*) growfromGM1->data;
    // nifti_image* growfromGM0 = nifti_image_read(fin, 1);
    // float* growfromGM0_data = (float*) growfromGM0->data;
    // growfromWM->dim[4] = 2;
    // nifti_update_dims_from_array(growfromGM);  // changing according sizes nt etc.

    nifti_image* GMkoordx1 = nifti_copy_nim_info(nim_input);
    nifti_image* GMkoordy1 = nifti_copy_nim_info(nim_input);
    nifti_image* GMkoordz1 = nifti_copy_nim_info(nim_input);
    nifti_image* GMkoordx2 = nifti_copy_nim_info(nim_input);
    nifti_image* GMkoordy2 = nifti_copy_nim_info(nim_input);
    nifti_image* GMkoordz2 = nifti_copy_nim_info(nim_input);

    GMkoordx1->datatype = NIFTI_TYPE_INT32;
    GMkoordy1->datatype = NIFTI_TYPE_INT32;
    GMkoordz1->datatype = NIFTI_TYPE_INT32;
    GMkoordx2->datatype = NIFTI_TYPE_INT32;
    GMkoordy2->datatype = NIFTI_TYPE_INT32;
    GMkoordz2->datatype = NIFTI_TYPE_INT32;

    GMkoordx1->nbyper = sizeof(int);
    GMkoordy1->nbyper = sizeof(int);
    GMkoordz1->nbyper = sizeof(int);
    GMkoordx2->nbyper = sizeof(int);
    GMkoordy2->nbyper = sizeof(int);
    GMkoordz2->nbyper = sizeof(int);

    GMkoordx1->data = calloc(GMkoordx1->nvox, GMkoordx1->nbyper);
    GMkoordy1->data = calloc(GMkoordy1->nvox, GMkoordy1->nbyper);
    GMkoordz1->data = calloc(GMkoordz1->nvox, GMkoordz1->nbyper);
    GMkoordx2->data = calloc(GMkoordx2->nvox, GMkoordx2->nbyper);
    GMkoordy2->data = calloc(GMkoordy2->nvox, GMkoordy2->nbyper);
    GMkoordz2->data = calloc(GMkoordz2->nvox, GMkoordz2->nbyper);

    int* GMkoordx1_data = (int*) GMkoordx1->data;
    int* GMkoordy1_data = (int*) GMkoordy1->data;
    int* GMkoordz1_data = (int*) GMkoordz1->data;
    int* GMkoordx2_data = (int*) GMkoordx2->data;
    int* GMkoordy2_data = (int*) GMkoordy2->data;
    int* GMkoordz2_data = (int*) GMkoordz2->data;

    nifti_image* equi_dist_layers  = nifti_copy_nim_info(nim_input);
    equi_dist_layers->datatype = NIFTI_TYPE_INT32;
    equi_dist_layers->nbyper = sizeof(int);
    equi_dist_layers->data = calloc(equi_dist_layers->nvox, equi_dist_layers->nbyper);
    int* equi_dist_layers_data = (int*) equi_dist_layers->data;

    // nifti_image* equi_dist_layers = nifti_image_read(fin, 1);
    // short  *equi_dist_layers_data = (short *) equi_dist_layers->data;
    // equi_dist_layers->dim[4] = 1;
    // nifti_update_dims_from_array(equi_dist_layers);  // changing according sizes nt etc.

    // Coordinates
    float x1g = 0.;
    float y1g = 0.;
    float z1g = 0.;

    float x2g = 0.;
    float y2g = 0.;
    float z2g = 0.;

    float x3g = 0.;
    float y3g = 0.;
    float z3g = 0.;

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

    int number_of_layers = 20;

    cout << "  Start growing from WM..." << endl;

    // Setting zero
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                *(growfromWM0_data + VOXEL_ID_3D) = 0;
                *(growfromWM1_data + VOXEL_ID_3D) = 0;

                *(growfromGM0_data + VOXEL_ID_3D) = 0;
                *(growfromGM1_data + VOXEL_ID_3D) = 0;
            }
        }
    }

    /////////////////////////
    // Closing 3D surfaces //
    /////////////////////////

    // cout << "  Closing surfaces..." << endl;
    //
    // int isWMb, isCSFb, isb;
    //
    // for (int iterations = 0; iterations< 100; ++iterations) {
    //     for (int iz=0; iz<sizeSlice; ++iz) {
    //         for (int iy=0; iy<sizePhase; ++iy) {
    //             for (int ix=0; ix<sizeRead-0; ++ix) {
    //                 if (*(nim_input_data + VOXEL_ID_3D) == 0) {
    //                     isWMb=0;
    //                     isCSFb=0;
    //                     isb=0;
    //                     for (int iz_i=max(0,iz-2); iz_i<min(sizeSlice,iz+2); ++iz_i) {
    //                         for (int iy_i=max(0,iy-2); iy_i<minxy * izn(sizePhase,iy+2); ++iy_i) {
    //                             for (int ix_i=max(0,ix-2); ix_i<min(sizeRead,ix+2); ++ix_i) {
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
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                if (*(nim_input_data + VOXEL_ID_3D) == 2) {
                    *(growfromWM0_data + VOXEL_ID_3D) = 1.;
                    *(WMkoordx1_data + VOXEL_ID_3D) = ix;
                    *(WMkoordy1_data + VOXEL_ID_3D) = iy;
                    *(WMkoordz1_data + VOXEL_ID_3D) = iz;
                }
            }
        }
    }

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead - 0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (*(nim_input_data + VOXEL_ID_3D) == 3 && *(growfromWM0_data + VOXEL_ID_3D) == 0) {
                        // cout << "  True   " << *(growfromWM0_data + nxy*islice + nx * ix + iy)<< endl;
                        for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc, sizeRead); ++ix_i) {
                                for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc, sizeRead); ++iz_i) {
                                    if (*(growfromWM0_data + nxy * iz_i + nx * ix_i + iy_i) == (float)grow_i) {
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
                            *(growfromWM0_data + VOXEL_ID_3D) = (float)grow_i+1;
                            *(WMkoordx1_data + VOXEL_ID_3D) = *(WMkoordx1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(WMkoordy1_data + VOXEL_ID_3D) = *(WMkoordy1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(WMkoordz1_data + VOXEL_ID_3D) = *(WMkoordz1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);

                        }
                        //cout << " ix   "  << ix << " iy   "  << iy  << "    " << *(WMkoordx1_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                    }
                }
            }
        }
    }

    ///////////////////
    // Grow from CSF //
    ///////////////////
    cout << "  Start growing from CSF ..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
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
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0.;
                    y1g = 0;
                    z1g = 0;
                    if (*(nim_input_data + VOXEL_ID_3D) == 3 && *(growfromGM0_data + VOXEL_ID_3D) == 0) {
                        for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc, sizeRead); ++iz_i) {
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc, sizePhase); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc, sizeRead); ++ix_i) {
                                    if (*(growfromGM0_data + nxy * iz_i + nx * ix_i + iy_i) == (float)grow_i) {
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
                            *(GMkoordx1_data + VOXEL_ID_3D) = *(GMkoordx1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(GMkoordy1_data + VOXEL_ID_3D) = *(GMkoordy1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(GMkoordz1_data + VOXEL_ID_3D) = *(GMkoordz1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
    }

    ////////////////////////////////////////////////////
    // Wabble across neigbouring voxles of closest WM //
    // to account for Pytagoras errors /////////////////
    ////////////////////////////////////////////////////
    cout << "  Correct for Pytagoras error..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                *(growfromWM1_data + VOXEL_ID_3D) = *(growfromWM0_data + VOXEL_ID_3D);
                *(WMkoordx2_data + VOXEL_ID_3D) = *(WMkoordx1_data + VOXEL_ID_3D);
                *(WMkoordy2_data + VOXEL_ID_3D) = *(WMkoordy1_data + VOXEL_ID_3D);
                *(WMkoordz2_data + VOXEL_ID_3D) = *(WMkoordz1_data + VOXEL_ID_3D);
            }
        }
    }

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    if (*(WMkoordy1_data + VOXEL_ID_3D) != 0) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        z1g = 0;

                        for (int iy_i = max(0, (*(WMkoordy2_data + VOXEL_ID_3D)) - grow_vinc); iy_i < min((*(WMkoordy2_data + VOXEL_ID_3D)) + grow_vinc, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, (*(WMkoordx2_data + VOXEL_ID_3D)) - grow_vinc); ix_i < min((*(WMkoordx2_data + VOXEL_ID_3D)) + grow_vinc, sizeRead); ++ix_i) {
                                for (int iz_i = max(0, (*(WMkoordz2_data + VOXEL_ID_3D)) - grow_vinc); iz_i < min((*(WMkoordz2_data + VOXEL_ID_3D)) + grow_vinc, sizeSlice); ++iz_i) {
                                    if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 2) {
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
                        *(WMkoordx2_data + VOXEL_ID_3D) = *(WMkoordx2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        *(WMkoordy2_data + VOXEL_ID_3D) = *(WMkoordy2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        *(WMkoordz2_data + VOXEL_ID_3D) = *(WMkoordz2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                    }
                }
            }
        }
    }
    cout << "  Running until stage 1..." << endl;

    //////////////////////////////////////////////////////
    // Wabble accross neigbouring voexles of closest GM //
    // to account for Pytagoras errors ///////////////////
    //////////////////////////////////////////////////////
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(growfromGM1_data + VOXEL_ID_3D) = *(growfromGM0_data + VOXEL_ID_3D);
                *(GMkoordx2_data + VOXEL_ID_3D) = *(GMkoordx1_data + VOXEL_ID_3D);
                *(GMkoordy2_data + VOXEL_ID_3D) = *(GMkoordy1_data + VOXEL_ID_3D);
                *(GMkoordz2_data + VOXEL_ID_3D) = *(GMkoordz1_data + VOXEL_ID_3D);
            }
        }
    }
    cout << "  Running until stage 2..." << endl;

    for (int grow_i = 1; grow_i < vinc; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix=0; ix < sizeRead - 0; ++ix) {
                    if (*(GMkoordy1_data + VOXEL_ID_3D) != 0) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        z1g = 0;
                        for (int iy_i = max(0, (*(GMkoordy2_data + VOXEL_ID_3D)) - grow_vinc); iy_i < min((*(GMkoordy2_data + VOXEL_ID_3D)) + grow_vinc, sizePhase); ++iy_i) {
                            for (int ix_i=max(0, (*(GMkoordx2_data + VOXEL_ID_3D)) - grow_vinc); ix_i < min((*(GMkoordx2_data + VOXEL_ID_3D)) + grow_vinc, sizeRead); ++ix_i) {
                                for (int iz_i=max(0, (*(GMkoordz2_data + VOXEL_ID_3D)) - grow_vinc); iz_i < min((*(GMkoordz2_data + VOXEL_ID_3D)) + grow_vinc, sizeRead); ++iz_i) {
                                    if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i)  == 1) {
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
                        *(GMkoordx2_data + VOXEL_ID_3D) = *(GMkoordx2_data + nxy*(int)z1g + nx * (int)x1g + (int)y1g);
                        *(GMkoordy2_data + VOXEL_ID_3D) = *(GMkoordy2_data + nxy*(int)z1g + nx * (int)x1g + (int)y1g);
                        *(GMkoordz2_data + VOXEL_ID_3D) = *(GMkoordz2_data + nxy*(int)z1g + nx * (int)x1g + (int)y1g);
                    }
                }
            }
        }
    }
    cout << "  Running until stage 3..." << endl;

    int GMK2_i, GMKz2_i, GMK3_i, WMK2_i, WMKz2_i, WMK3_i;
    float GMK2_f, GMKz2_f, GMK3_f, WMK2_f, WMKz2_f, WMK3_f, ix_f, iy_f, iz_f;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
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
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
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
    return sqrt((x1 - x2) * (x1 - x2) * dX * dX + (y1 - y2) * (y1 - y2)
                * dY * dY + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float angle(float a, float b, float c) {
    if (a * a + b * b - c * c <= 0) {
        return 3.141592;
    } else {
        return acos((a * a + b * b - c * c) / (2. * a * b));
    }
}
