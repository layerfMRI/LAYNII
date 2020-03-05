
// TODO(@Faruk): Test requires `layers.nii` I need to get it from Renzo.


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_3DCOLUMNS : Calculates cortical distances (columnar structures) \n"
    "               based on the gray matter (GM) geometry.\n"
    "\n"
    "Usage:\n"
    "    LN_3DCOLUMNS -layer_file layers.nii -landmarks landmarks.nii \n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -layer_file   : Nifti (.nii) file containing layer or column masks \n"
    "    -landmarks    : Nifti (.nii) file with landmarks 1, 2, 3 (1 is in \n"
    "                    the center 2 and 3 are the borders. Landmarks \n"
    "                    should be at least 4 voxels thick.\n"
    "    -twodim       : (Optional) Run in 2D only.\n"
    "    -mask         : (Optional) Mask activity outside of layers.\n"
    "    -vinc         : Number of columns.\n"
    "    -jiajiaoption : Include cerebrospinal fluid (CSF). Only do this \n"
    "                    if two sides of the sulcus are not touching. \n"
    "\n"
    "Notes:\n"
    "     - Layer nifti and landmarks nifti should have the same dimensions \n"
    "     - This program now supports INT16, INT32 and FLOAT23 \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * layer_filename = NULL, *fout = NULL, *landmarks_filename = NULL;
    int ac, twodim = 0, do_masking = 0, jiajiavinc_max = 45, jiajiaoption = 0;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
    }
    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layer_file\n");
                return 1;
            }
            layer_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-landmarks")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            landmarks_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << "Smoothing only in 2D." << endl;
        } else if (!strcmp(argv[ac], "-jiajiaoption")) {
            jiajiaoption = 1;
            cout << "Do not remove CSF." << endl;
        } else if (!strcmp(argv[ac], "-vinc")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -vinc\n");
                return 1;
            }
            jiajiavinc_max = atof(argv[ac]);  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-mask")) {
            do_masking = 1;
            cout << "Set every voxel to zero outside the layers." << endl;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!landmarks_filename) {
        fprintf(stderr, "** missing option '-landmarks'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_landmarks_r = nifti_image_read(landmarks_filename, 1);
    if (!nim_landmarks_r) {
        fprintf(stderr, "** failed to read landmarks NIfTI from '%s'\n", landmarks_filename);
        return 2;
    }
    if (!layer_filename) { fprintf(stderr, "** missing option '-layer_file'\n");  return 1; }
    // Read input dataset, including data
    nifti_image * nim_layers_r = nifti_image_read(layer_filename, 1);
    if (!nim_layers_r) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", layer_filename);
        return 2;
    }

    log_welcome("LN_3DCOLUMNS");
    log_nifti_descriptives(nim_layers_r);
    log_nifti_descriptives(nim_landmarks_r);

    // Get dimensions of input
    int sizeSlice = nim_layers_r->nz;
    int sizePhase = nim_layers_r->nx;
    int sizeRead = nim_layers_r->ny;
    int nrep = nim_layers_r->nt;
    int nx = nim_layers_r->nx;
    int nxy = nim_layers_r->nx * nim_layers_r->ny;
    int nxyz = nim_layers_r->nx * nim_layers_r->ny * nim_layers_r->nz;
    float dX = nim_layers_r->pixdim[1];
    float dY = nim_layers_r->pixdim[2];
    float dZ = nim_layers_r->pixdim[3];

    if (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_layers = nifti_copy_nim_info(nim_layers_r);
    nim_layers->datatype = NIFTI_TYPE_FLOAT32;
    nim_layers->nbyper = sizeof(float);
    nim_layers->data = calloc(nim_layers->nvox, nim_layers->nbyper);
    float *nim_layers_data = (float *) nim_layers->data;

    nifti_image * nim_landmarks = nifti_copy_nim_info(nim_layers_r);
    nim_landmarks->datatype = NIFTI_TYPE_INT32;
    nim_landmarks->nbyper = sizeof(int);
    nim_landmarks->data = calloc(nim_landmarks->nvox, nim_landmarks->nbyper);
    int *nim_landmarks_data = (int *) nim_landmarks->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_landmarks_r->datatype == NIFTI_TYPE_FLOAT32) {
        float  *nim_landmarks_r_data = (float *) nim_landmarks_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_landmarks_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_landmarks_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_landmarks_r->datatype == NIFTI_TYPE_INT16) {
        short  *nim_landmarks_r_data = (short *) nim_landmarks_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_landmarks_data + nxyz *it + nxy*islice + nx*ix + iy) = (float) (*(nim_landmarks_r_data + nxyz *it + nxy*islice + nx*ix + iy));
                    }
                }
            }
        }
    }
    if (nim_landmarks_r->datatype == NIFTI_TYPE_INT32) {
        int *nim_landmarks_r_data = (int *) nim_landmarks_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_landmarks_data + nxyz *it + nxy * islice + nx * ix + iy) = (float) (*(nim_landmarks_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_FLOAT32) {
        float  *nim_layers_r_data = (float *) nim_layers_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_INT16) {
        short  *nim_layers_r_data = (short *) nim_layers_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_INT32) {
        int *nim_layers_r_data = (int *) nim_layers_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    //////////////////////////////
    // Finding number of layers //
    //////////////////////////////
    int layernumber = 0;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > layernumber) layernumber = *(nim_layers_data + nxy * iz + nx * ix + iy);
            }
        }
    }
    cout << "  There are  " << layernumber<< " layers." << endl;

    ////////////////////////////////
    // Allocating necessary files //
    ////////////////////////////////
    nifti_image * Grow_x = nifti_copy_nim_info(nim_layers);
    nifti_image * Grow_y = nifti_copy_nim_info(nim_layers);
    nifti_image * Grow_z = nifti_copy_nim_info(nim_layers);

    Grow_x->datatype = NIFTI_TYPE_INT16;
    Grow_y->datatype = NIFTI_TYPE_INT16;
    Grow_x->datatype = NIFTI_TYPE_INT16;

    Grow_x->nbyper = sizeof(short);
    Grow_y->nbyper = sizeof(short);
    Grow_z->nbyper = sizeof(short);

    Grow_x->data = calloc(Grow_x->nvox, Grow_x->nbyper);
    Grow_y->data = calloc(Grow_y->nvox, Grow_y->nbyper);
    Grow_z->data = calloc(Grow_z->nvox, Grow_z->nbyper);

    short *Grow_x_data = (short *) Grow_x->data;
    short *Grow_y_data = (short *) Grow_y->data;
    short *Grow_z_data = (short *) Grow_z->data;

    // cout << "  Not Segmentation fault until here 272" << endl;

    nifti_image * growfromCenter = nifti_copy_nim_info(nim_layers);
    nifti_image * growfromCenter_thick = nifti_copy_nim_info(nim_layers);
    nifti_image * growfromLeft = nifti_copy_nim_info(nim_layers);
    nifti_image * growfromRight = nifti_copy_nim_info(nim_layers);

    growfromCenter->datatype = NIFTI_TYPE_INT16;
    growfromCenter_thick->datatype = NIFTI_TYPE_INT16;
    growfromLeft->datatype = NIFTI_TYPE_INT16;
    growfromRight->datatype = NIFTI_TYPE_INT16;

    growfromCenter->nbyper = sizeof(short);
    growfromCenter_thick->nbyper = sizeof(short);
    growfromLeft->nbyper = sizeof(short);
    growfromRight->nbyper = sizeof(short);

    growfromCenter->data = calloc(growfromCenter->nvox, growfromCenter->nbyper);
    growfromCenter_thick->data = calloc(growfromCenter_thick->nvox, growfromCenter_thick->nbyper);
    growfromLeft->data = calloc(growfromLeft->nvox, growfromLeft->nbyper);
    growfromRight->data = calloc(growfromRight->nvox, growfromRight->nbyper);

    short  *growfromCenter_data = (short *) growfromCenter->data;
    short  *growfromCenter_thick_data = (short *) growfromCenter_thick->data;
    short  *growfromLeft_data = (short *) growfromLeft->data;
    short  *growfromRight_data = (short *) growfromRight->data;

    // cout << "  Not Segmentation fault until here 292" << endl;

    nifti_image * lateralCoord = nifti_copy_nim_info(nim_layers);
    nifti_image * inferioCoord = nifti_copy_nim_info(nim_layers);
    lateralCoord->datatype = NIFTI_TYPE_INT16;
    inferioCoord->datatype = NIFTI_TYPE_INT16;
    lateralCoord->nbyper = sizeof(short);
    inferioCoord->nbyper = sizeof(short);
    lateralCoord->data = calloc(lateralCoord->nvox, lateralCoord->nbyper);
    inferioCoord->data = calloc(inferioCoord->nvox, inferioCoord->nbyper);
    short *lateralCoord_data = (short *) lateralCoord->data;
    short *inferioCoord_data = (short *) inferioCoord->data;

    // cout << "  Not Segmentation fault until here 306" << endl;

    ///////////////////////////////
    // Prepare growing variables //
    ///////////////////////////////
    float x1g = 0.;
    float y1g = 0.;
    float z1g = 0.;

    float dist_min2 = 0.;
    float dist_i = 0.;
    float dist_p1 = 0.;

    int grow_vinc = 3;
    int grow_vinc_area = 1;
    int vinc_max = 250;
    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);

    /////////////////////////
    // Growing from Center //
    /////////////////////////
    cout << "  Growing from center " << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                // Defining seed at center landmark
                if (*(nim_landmarks_data + nxy * iz + nx * ix + iy) == 1 && abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber / 2)) < 2) {
                    *(growfromCenter_data + nxy * iz + nx * ix + iy) = 1.;
                    *(Grow_x_data + nxy * iz + nx * ix + iy) = ix;
                    *(Grow_y_data + nxy * iz + nx * ix + iy) = iy;
                    *(Grow_z_data + nxy * iz + nx * ix + iy) = iz;
                }
            }
        }
    }
    for (int grow_i = 1; grow_i < vinc_max; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber/2)) < 2 && *(growfromCenter_data + nxy * iz + nx * ix + iy) == 0 && *(nim_landmarks_data + nxy * iz + nx * ix + iy) < 2) {
                        // NOTE: Only grow into areas that are GM and that have not been grown into, yet.
                        // And it should stop as soon as it hits tie border.
                        for (int iy_i = max(0, iy - grow_vinc_area); iy_i<min(iy + grow_vinc_area + 1, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, ix - grow_vinc_area); ix_i < min(ix + grow_vinc_area + 1, sizeRead); ++ix_i) {
                                for (int iz_i = max(0, iz - grow_vinc_area); iz_i < min(iz + grow_vinc_area + 1, sizeSlice); ++iz_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);

                                    if (*(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i) == grow_i && *(nim_landmarks_data + nxy * iz + nx * ix + iy) < 2) {
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
                        if (dist_min2 < 1.7) {  // TODO(@Renzo): ???? I DONT REMEMBER WHY I NEED THIS ????
                            // distDebug(0, islice, iy, ix) = dist_min2;
                            *(growfromCenter_data + nxy * iz + nx * ix + iy) = grow_i+1;
                            *(Grow_x_data + nxy * iz + nx * ix + iy) = *(Grow_x_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_y_data + nxy * iz + nx * ix + iy) = *(Grow_y_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_z_data + nxy * iz + nx * ix + iy) = *(Grow_z_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                        //cout << " ix   " << ix << " iy   " << iy << "    " << *(WMkoord0_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                    }
                }
            }
        }
    }

    const char *fout_5 = "finding_leaks.nii";
    if (nifti_set_filenames(growfromCenter, fout_5, 1, 1)) return 1;
    nifti_image_write(growfromCenter);

    ///////////////////////
    // Growing from left //
    ///////////////////////
    cout << "  Growing from left..." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                *(Grow_x_data + nxy * iz + nx * ix + iy) = 0;
                *(Grow_y_data + nxy * iz + nx * ix + iy) = 0;
                *(Grow_z_data + nxy * iz + nx * ix + iy) = 0;
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                // Defining seed at center landmark
                if (*(nim_landmarks_data + nxy * iz + nx * ix + iy) == 2 && abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber / 2)) < 2) {
                    *(growfromLeft_data + nxy * iz + nx * ix + iy) = 1.;
                    *(Grow_x_data + nxy * iz + nx * ix + iy) = ix;
                    *(Grow_y_data + nxy * iz + nx * ix + iy) = iy;
                    *(Grow_z_data + nxy * iz + nx * ix + iy) = iz;
                }
            }
        }
    }
    for (int grow_i = 1; grow_i < vinc_max; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber / 2)) < 2 && *(growfromLeft_data + nxy * iz + nx * ix + iy) == 0  && *(growfromCenter_data + nxy * iz + nx * ix + iy) != 0) {
                        // NOTE: Only grow into areas that are GM and that have not been grown into, yet ...
                        // And it should stop as soon as it hits the border.
                        for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, sizeRead); ++ix_i) {
                                for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc + 1, sizeSlice); ++iz_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);

                                    if (*(growfromLeft_data + nxy * iz_i + nx * ix_i + iy_i) == grow_i) {
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
                        if (dist_min2 < 1.7) {  // TODO(@Renzo): ???? I DONT REMEMBER WHY I NEED THIS ????
                            // distDebug(0, islice, iy, ix) = dist_min2;
                            *(growfromLeft_data + nxy * iz + nx * ix + iy) = grow_i + 1;
                            *(Grow_x_data + nxy * iz + nx * ix + iy) = *(Grow_x_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_y_data + nxy * iz + nx * ix + iy) = *(Grow_y_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_z_data + nxy * iz + nx * ix + iy) = *(Grow_z_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                        //cout << " ix   " << ix << " iy   " << iy << "    " << *(WMkoord0_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                    }
                }
            }
        }
    }

    ////////////////////////
    // Growing from Right //
    ////////////////////////
    cout << "  Growing from right..." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(Grow_x_data + nxy*iz + nx*ix + iy) = 0;
                *(Grow_y_data + nxy*iz + nx*ix + iy) = 0;
                *(Grow_z_data + nxy*iz + nx*ix + iy) = 0;
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                // Defining seed at center landmark
                if (*(nim_landmarks_data + nxy * iz + nx * ix + iy) == 3 && abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber / 2)) < 2) {
                    *(growfromRight_data + nxy * iz + nx * ix + iy) = 1.;
                    *(Grow_x_data + nxy * iz + nx * ix + iy) = ix;
                    *(Grow_y_data + nxy * iz + nx * ix + iy) = iy;
                    *(Grow_z_data + nxy * iz + nx * ix + iy) = iz;
                }
            }
        }
    }
    for (int grow_i = 1; grow_i < vinc_max; grow_i++) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead - 0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber/2)) < 2 && *(growfromRight_data + nxy * iz + nx * ix + iy) == 0  && *(growfromCenter_data + nxy * iz + nx * ix + iy) != 0) {
                        // NOTE: Only grow into areas that are GM and that have not been gown into, yet...
                        // And it should stop as soon as it hits tie border
                        for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, sizeRead); ++ix_i) {
                                for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc + 1, sizeSlice); ++iz_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);

                                    if (*(growfromRight_data + nxy * iz_i + nx * ix_i + iy_i) == grow_i) {
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
                        if (dist_min2 < 1.7) {  // TODO(@Renzo): ???? I DONT REMEMBER WHY I NEED THIS ????
                            // distDebug(0, islice, iy, ix) = dist_min2;
                            *(growfromRight_data + nxy * iz + nx * ix + iy) = grow_i+1;
                            *(Grow_x_data + nxy * iz + nx * ix + iy) = *(Grow_x_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_y_data + nxy * iz + nx * ix + iy) = *(Grow_y_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_z_data + nxy * iz + nx * ix + iy) = *(Grow_z_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                        // cout << " ix   " << ix << " iy   " << iy << "    " << *(WMkoord0_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                    }
                }
            }
        }
    }

    const char  *fout_11 = "grow_from_right.nii";
    if (nifti_set_filenames(growfromRight, fout_11, 1, 1)) return 1;
    nifti_image_write(growfromRight);

    //////////////////////////////////////
    // Get normalized coordinate system //
    //////////////////////////////////////
    cout << "  Getting normalized coordinate system..." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                if (*(growfromCenter_data + nxy * iz + nx * ix + iy) > 0) {
                    *(lateralCoord_data + nxy * iz + nx * ix + iy) = (int) (200 *(float) *(growfromRight_data + nxy * iz + nx * ix + iy) / ((float) *(growfromRight_data + nxy * iz + nx * ix + iy) + (float) *(growfromLeft_data + nxy * iz + nx * ix + iy)));
                    // *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) = 1;
                }
            }
        }
    }

    const char  *fout_9 = "lateral_cord.nii";
    if (nifti_set_filenames(lateralCoord, fout_9, 1, 1)) return 1;
    nifti_image_write(lateralCoord);

    //////////////////////////
    // Sparse visualisation //
    //////////////////////////
    // cout << " sparse visualisation " << endl;
    // cout << " visualisation " << endl;
    // for (int iz = 0; iz < sizeSlice; ++iz) {
    //     for (int iy = 0; iy < sizePhase; ++iy) {
    //         for (int ix = 0; ix < sizeRead - 0; ++ix) {
    //             if (iy % 4 == 0 || iz % 4 == 0) {
    //                 *(lateralCoord_data + nxy * iz + nx * ix + iy) = *(lateralCoord_data + nxy * iz + nx * ix + iy);
    //             } else {
    //                 *(lateralCoord_data + nxy * iz + nx * ix + iy) = 0;
    //             }
    //         }
    //     }
    // }
    // cout << "  Not Segmentation fault until here 381" << endl;

    ////////////////////
    // Smooth columns //
    ////////////////////
    // cout << "  Smooth columns " << endl;

    float gaus(float distance, float sigma);
    cout << "  Smoothing in middle layer..." << endl;

    nifti_image * smoothed = nifti_copy_nim_info(nim_layers);
    nifti_image * gausweight = nifti_copy_nim_info(nim_layers);
    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
    float *smoothed_data = (float *) smoothed->data;
    float *gausweight_data = (float *) gausweight->data;

    // float kernal_size = 10; // Corresponds to one voxel size.
    int FWHM_val = 1;
    int vinc_sm = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
    dist_i = 0.;
    cout << "  vinc_sm " << vinc_sm<< endl;
    cout << "  FWHM_val " << FWHM_val<< endl;
    // cout << " DEBUG " << dist(1., 1., 1., 1., 2., 1., dX, dY, dZ) << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;

                if (*(lateralCoord_data + nxy * iz + nx * ix + iy)  > 0) {
                    for (int iz_i = max(0, iz - vinc_sm); iz_i < min(iz + vinc_sm + 1, sizeSlice); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc_sm); iy_i < min(iy + vinc_sm + 1, sizePhase); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc_sm); ix_i<min(ix + vinc_sm + 1, sizeRead); ++ix_i) {
                                if (*(lateralCoord_data + nxy * iz_i + nx * ix_i + iy_i)  > 0) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    // cout << "debug 4 " << gaus(dist_i, FWHM_val) << endl;
                                    // cout << "debug 5 " << dist_i << endl;
                                    // if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) cout << "debug  4b " << endl;
                                    // dummy = *(layer_data + nxy * iz_i + nx * ix_i + iy_i);
                                    *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(lateralCoord_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                            }
                        }
                    }
                    if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) {
                        *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) / *(gausweight_data + nxy * iz + nx * ix + iy);
                    }
                }
                //if (*(nim_layers_r_data + nxy * iz + nx * ix + iy) <= 0) *(smoothed_data + nxy*iz + nx * ix + iy) = *(lateralCoord_data + nxy*iz + nx*ix + iy);
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(growfromCenter_data + nxy * iz + nx * ix + iy) > 0) {
                    *(lateralCoord_data + nxy * iz + nx * ix + iy) = (int) *(smoothed_data + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    // const char *fout_10 = "lateral_cord_smoothed.nii";
    // if (nifti_set_filenames(lateralCoord, fout_10, 1, 1)){
    //     return 1;
    // }
    // nifti_image_write(lateralCoord);

    /////////////////////////////////////
    // Extending columns across layers //
    /////////////////////////////////////
    cout << "  Extending columns across layers..." << endl;

    nifti_image * hairy_brain = nifti_copy_nim_info(nim_layers);
    hairy_brain->datatype = NIFTI_TYPE_INT16;
    hairy_brain->nbyper = sizeof(short);
    hairy_brain->data = calloc(hairy_brain->nvox, hairy_brain->nbyper);
    short *hairy_brain_data = (short *) hairy_brain->data;
    nifti_image * hairy_brain_dist = nifti_copy_nim_info(nim_layers);
    hairy_brain_dist->datatype = NIFTI_TYPE_FLOAT32;
    hairy_brain_dist->nbyper = sizeof(float);
    hairy_brain_dist->data = calloc(hairy_brain_dist->nvox, hairy_brain_dist->nbyper);
    float  *hairy_brain_dist_data = (float *) hairy_brain_dist->data;

    dist_min2 = 10000.;  // This is an upper limit of the cortical thickness
    int vinc_thickness = 13;  // Closest middle layer area algorithm looks for
    int cloasest_coord = 0.;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                *(hairy_brain_dist_data + nxy * iz + nx * ix + iy) = dist_min2;
            }
        }
    }

    // *(nim_layers_data + nxy * iz + nx * ix + iy) > 0

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                dist_min2 = 10000.;
                if (*(lateralCoord_data + nxy*iz + nx*ix + iy) > 0) {
                    // for (int iz_i = max(0, iz-vinc_thickness + 1); iz_i < min(iz + vinc_thickness + 1, sizeSlice); ++iz_i) {
                    int iz_i = iz;
                    for (int iy_i = max(0, iy-vinc_thickness+1); iy_i<min(iy+vinc_thickness+1, sizePhase); ++iy_i) {
                        for (int ix_i = max(0, ix-vinc_thickness+1); ix_i<min(ix+vinc_thickness+1, sizeRead); ++ix_i) {
                            if (*(lateralCoord_data + nxy*iz_i + nx*ix_i + iy_i) == 0  && *(nim_layers_data + nxy*iz + nx*ix + iy) > 0 && dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ) < vinc_thickness) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);

                                if (dist_i < *(hairy_brain_dist_data + nxy * iz_i + nx * ix_i + iy_i) && *(nim_layers_data + nxy * iz + nx * ix + iy) > 0) {
                                    *(hairy_brain_dist_data + nxy * iz_i + nx * ix_i + iy_i) = dist_i;
                                    *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = *(lateralCoord_data + nxy * iz + nx * ix + iy);
                                }
                            }
                            if (*(lateralCoord_data + nxy * iz_i + nx * ix_i + iy_i) != 0) {
                                *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = *(lateralCoord_data + nxy * iz_i + nx * ix_i + iy_i);
                            }
                        }
                    }
                    // }
                    // *(hairy_brain_data + nxy * iz + nx * ix + iy) = cloasest_coord;
                }
            }
        }
    }

    // const char *fout_7 = "thick_ribbon.nii";
    // if (nifti_set_filenames(hairy_brain, fout_7, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(hairy_brain);

    ////////////////////
    // Smooth columns //
    ////////////////////
    cout << "  Smoothing hairy brain again... " << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                *(smoothed_data + nxy * iz + nx *ix + iy) = 0;
            }
        }
    }

    FWHM_val = 1;
    vinc_sm = max(1., 2. * FWHM_val / dX);
    cout << "  vinc_sm " << vinc_sm << endl;
    cout << "  FWHM_val " << FWHM_val << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead -0; ++ix) {
                *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;

                if (*(hairy_brain_data + nxy * iz + nx * ix + iy) > 0) {
                    int layernumber_i = *(nim_layers_data + nxy * iz + nx * ix + iy);

                    //for (int iz_i = max(0, iz - vinc_sm); iz_i < min(iz + vinc_sm + 1, sizeRead); ++iz_i) {
                    int iz_i = iz;
                    for (int iy_i = max(0, iy - vinc_sm); iy_i < min(iy + vinc_sm + 1, sizePhase); ++iy_i) {
                        for (int ix_i = max(0, ix - vinc_sm); ix_i < min(ix + vinc_sm + 1, sizeRead); ++ix_i) {
                            if (abs((int) *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) - layernumber_i) < 2 && *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) > 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                // cout << "debug  4 " << gaus(dist_i, FWHM_val) << endl;
                                // cout << "debug  5 " << dist_i << endl;
                                // if (*(nim_input_data + nxy*iz + nx*ix + iy) == 3) cout << "debug  4b " << endl;
                                // dummy = *(layer_data + nxy*iz_i + nx*ix_i + iy_i);
                                *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                            }
                        }
                        // }
                    }
                    if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) {
                        *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) / *(gausweight_data + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(hairy_brain_data + nxy * iz + nx * ix + iy) > 0) {
                    *(hairy_brain_data + nxy * iz + nx * ix + iy) = (int) *(smoothed_data + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    // const char *fout_8 = "thick_ribbon_smoothed.nii";
    // if (nifti_set_filenames(hairy_brain, fout_8, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(hairy_brain);

    //////////////////////////
    // Sparse visualisation //
    //////////////////////////
    cout << "  Visualisation..." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                if ((int) *(smoothed_data + nxy * iz + nx * ix + iy) % 10 == 0 && iz % 10 == 0) {
                    *(growfromRight_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx *ix + iy);
                }
                else *(growfromRight_data + nxy * iz + nx * ix + iy) = 0;
            }
        }
    }
    // cout << "  Not Segmentation fault until here 381" << endl;

    //////////////////////////////////
    // Growing from as thick cortex //
    //////////////////////////////////
    cout << "  Growing from center with thick cortex..." << endl;
    cout << "  Jiajia Option is on." << endl;

    int Jiajia_otion = 1;
    int grow_vinc_thick = 1;
    int vinc_max_thick = 17;
    int grow_vinc_area_thick = 1;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                // Defining seed at center landmark
                if (*(growfromCenter_data + nxy * iz + nx * ix + iy) > 0) {
                    *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) = 1;
                }
            }
        }
    }
    if (jiajiaoption == 1)  {
        for (int grow_i = 1; grow_i < vinc_max_thick; grow_i++) {
            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead-0; ++ix) {
                        if ( *(nim_layers_data + nxy * iz + nx * ix + iy) > 0 && *(nim_layers_data + nxy * iz + nx * ix + iy) <= layernumber && *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) == grow_i && *(hairy_brain_data + nxy * iz + nx *ix + iy) > 0) {
                            // Note: Only grow into areas that are GM and that have not been grown into, yet...
                            // And it should stop as soon as it hits the border
                            int iz_i = iz;
                            // int stopme = 0;
                            for (int ix_i = max(0, ix - grow_vinc_area_thick); ix_i < min(ix + grow_vinc_area_thick + 1, sizeRead); ++ix_i) {
                                for (int iy_i = max(0, iy - grow_vinc_area_thick); iy_i < min(iy + grow_vinc_area_thick + 1, sizePhase); ++iy_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    // if (*(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) == layernumber) stopme = 1;
                                    if (dist_i <= (dY + dX) / 2. && *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) == 0 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) <= layernumber && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) > 0 && *(hairy_brain_data + nxy * iz + nx * ix + iy) > 0) {
                                        *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) = grow_i+1;
                                    }
                                }
                            }
                            // cout << " ix   " << ix << " iy   " << iy << "    " << *(WMkoord0_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                        }
                    }
                }
            }
        }
    }
    if (jiajiaoption != 1) {
        for (int grow_i = 1; grow_i < vinc_max_thick; grow_i++) {
            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy <sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead-0; ++ix) {
                        if ( *(nim_layers_data + nxy*iz + nx*ix + iy) > 0 && *(nim_layers_data + nxy*iz + nx*ix + iy) < layernumber && *(growfromCenter_thick_data + nxy*iz + nx*ix + iy) == grow_i && *(hairy_brain_data + nxy*iz + nx*ix + iy) > 0) {
                            // Note: Only grow into areas that are GM and that have not been gown into, yet...
                            // And it should stop as soon as it hits the border.
                            int iz_i = iz;
                            // int stopme = 0;
                            for (int ix_i = max(0, ix - grow_vinc_area_thick); ix_i < min(ix + grow_vinc_area_thick + 1, sizeRead); ++ix_i) {
                                for (int iy_i = max(0, iy - grow_vinc_area_thick); iy_i<min(iy + grow_vinc_area_thick + 1, sizePhase); ++iy_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    // if (*(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) == layernumber) stopme = 1;
                                    if (dist_i <= (dY + dX) / 2. && *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) == 0 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) < layernumber && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) > 0 && *(hairy_brain_data + nxy * iz + nx * ix + iy) > 0) {
                                        *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) = grow_i + 1;
                                    }
                                }
                            }
                            // cout << " ix   " << ix << " iy   " << iy << "    " << *(WMkoord0_data + nxy*islice + nx*(int)x1g + (int)y1g)<< endl;
                        }
                    }
                }
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                if ((*(nim_layers_data + nxy * iz + nx * ix + iy) == layernumber - 1 || *(nim_layers_data + nxy * iz + nx * ix + iy) == layernumber - 2) && *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) > 0) {
                    // only grow into areas that are GM and that have not been gown into, yet .... and it should stop as soon as it hits tie border
                    for (int iy_i = max(0, iy - grow_vinc_area_thick); iy_i < min(iy + grow_vinc_area_thick + 1, sizePhase); ++iy_i) {
                        for (int ix_i = max(0, ix - grow_vinc_area_thick); ix_i < min(ix + grow_vinc_area_thick + 1, sizeRead); ++ix_i) {
                            int iz_i = iz;
                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);

                            if (dist_i <= (dY + dX) / 2. && *(nim_landmarks_data + nxy * iz_i + nx * ix_i + iy_i) > 0) {
                                *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) = *(growfromCenter_thick_data + nxy * iz + nx * ix + iy);
                            }
                        }
                    }
                }
            }
        }
    }
    // const char  *fout_6 = "growfromCenter_thick.nii";
    // if (nifti_set_filenames(growfromCenter_thick, fout_6, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(growfromCenter_thick);

    /////////////////////////////
    // Clean up hairy brain /////
    /////////////////////////////
    cout << "  Cleaning up hairy brain..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                if (*(growfromCenter_thick_data + nxy * iz + nx * ix + iy) == 0 || *(hairy_brain_data + nxy * iz + nx * ix + iy) == 0) {
                    *(growfromRight_data + nxy*iz + nx*ix + iy) = 0;
                    *(hairy_brain_data + nxy*iz + nx*ix + iy) = 0;
                }
            }
        }
    }

    //////////////////////////////
    // Change number of columns //
    //////////////////////////////
    cout << "  Number of columns." << endl;
    int max_columns = 0;
    int min_columns = 100000000;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                if (*(hairy_brain_data + nxy * iz + nx * ix + iy) >  0) {
                    if ((int) *(hairy_brain_data + nxy * iz + nx *ix + iy) > max_columns) max_columns = (int) *(hairy_brain_data + nxy * iz + nx * ix + iy);
                    if ((int) *(hairy_brain_data + nxy * iz + nx *ix + iy) < min_columns) min_columns = (int) *(hairy_brain_data + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    cout << "  Max = " << max_columns << " | Min = " << min_columns << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead - 0; ++ix) {
                if (*(hairy_brain_data + nxy * iz + nx * ix + iy) >  0) {
                    *(hairy_brain_data + nxy * iz + nx * ix + iy) = (*(hairy_brain_data + nxy * iz + nx * ix + iy) - (short)min_columns) * (short)jiajiavinc_max / (short)(max_columns - min_columns);
                }
            }
        }
    }
    // cout << "  Not Segmentation fault until here 381" << endl;

    cout << "  Writing output... " << endl;

    // // Output file name
    // const char *fout_4 = "leaky_layers.nii";
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(leak_layer);

    // const char *fout_2 = "smoothed_coordinates.nii";
    // if (nifti_set_filenames(smoothed, fout_2, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(smoothed);

    // const char *fout_4 = "growfromLeft.nii";
    // if (nifti_set_filenames(growfromLeft, fout_4, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(growfromLeft);

    // const char *fout_5 = "finding_leaks.nii";
    // if (nifti_set_filenames(growfromCenter, fout_5, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(growfromCenter);

    // const char *fout_3 = "lateral_coord.nii";
    // if (nifti_set_filenames(lateralCoord, fout_3, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(lateralCoord);

    const char  *fout_1 = "column_coordinates.nii";
    if (nifti_set_filenames(hairy_brain, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(hairy_brain);

    cout << "  Finished." << endl;
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt((x1 - x2) * (x1 - x2) * dX * dX
                + (y1 - y2) * (y1 - y2) * dY * dY
                + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592)) * exp(-0.5 * distance * distance / (sigma * sigma));
}
