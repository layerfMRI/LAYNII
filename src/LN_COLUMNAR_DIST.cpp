
// TODO(Faruk): Requires layers.nii from Renzo for testing.


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_COLUMNAR_DIST : Calculates cortical distances (columnar structures) \n"
    "                   based on the gray matter geometry."
    "\n"
    "Usage:\n"
    "    LN_COLUMNAR_DIST -layer_file layers.nii -landmarks landmarks.nii \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -layer_file : Nifti (.nii) that contains layer or column masks. \n"
    "    -landmarks  : Nifti (.nii) with landmarks (use value 1 as origin). \n"
    "                  Landmarks should be at least 4 voxels thick. \n"
    "    -vinc       : Maximal length of cortical distance. Bigger values \n"
    "                  will take longer but span a larger cortical area. \n"
    "                  Default is 40. \n"
    "    -Ncolumns   : (Optional) For the number of columns. Smaller values \n"
    "                  will result in thick columns \n"
    "    -twodim     : (Optional) Run in 2D only. Though this will not \n"
    "                  really make it faster. \n"
    "    -verbose    : (Optional) to write out all the intermediate \n"
    "                  steps of the algorithm (e.g. for debugging) \n"
    "\n"
    "Notes:\n"
    "    - The layer nii file and the landmarks nii file should have the \n"
    "      same dimensions.\n"
    "    - This program now supports INT16, INT32 and FLOAT32. \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * layer_filename = NULL, * fout = NULL, * landmarks_filename = NULL;
    int ac, twodim = 0, do_masking = 0, vinc_max = 40, Ncolumns = 0, verbose = 0;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -layer_file\n");
                return 1;
            }
            layer_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-landmarks")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -input\n");
                return 1;
            }
            landmarks_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << " I will do smoothing only in 2D" << endl;
        } else if (!strcmp(argv[ac], " - vinc")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for  - vinc\n");
                return 1;
            }
            vinc_max = atof(argv[ac]);  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-Ncolumns")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -Ncolumns\n");
                return 1;
            }
            Ncolumns = atof(argv[ac]);  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-verbose")) {
            verbose = 1;
            cout << " I will give you everything I have, happy debugging" << endl;
        } else {
            fprintf(stderr, " ** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!landmarks_filename) {
        fprintf(stderr, " ** missing option '-landmarks'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_landmarks_r = nifti_image_read(landmarks_filename, 1);
    if (!nim_landmarks_r) {
        fprintf(stderr, " ** failed to read layer NIfTI image from '%s'\n", landmarks_filename);
        return 2;
    }

    if (!layer_filename) { fprintf(stderr, " ** missing option '-layer_file'\n");  return 1; }
    // read input dataset, including data
    nifti_image * nim_layers_r = nifti_image_read(layer_filename, 1);
    if (!nim_layers_r) {
        fprintf(stderr, " ** failed to read layer NIfTI image from '%s'\n", layer_filename);
        return 2;
    }

    log_welcome("LN_COLUMNAR_DIST");
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

    if  (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_layers = nifti_copy_nim_info(nim_layers_r);
    nim_layers->datatype = NIFTI_TYPE_FLOAT32;
    nim_layers->nbyper = sizeof(float);
    nim_layers->data = calloc(nim_layers->nvox, nim_layers->nbyper);
    float * nim_layers_data = (float *) nim_layers->data;

    nifti_image * nim_landmarks = nifti_copy_nim_info(nim_layers_r);
    nim_landmarks->datatype = NIFTI_TYPE_INT32;
    nim_landmarks->nbyper = sizeof(int);
    nim_landmarks->data = calloc(nim_landmarks->nvox, nim_landmarks->nbyper);
    int * nim_landmarks_data = (int *) nim_landmarks->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_landmarks_r->datatype == NIFTI_TYPE_FLOAT32) {
        float * nim_landmarks_r_data = (float *) nim_landmarks_r->data;
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
        short * nim_landmarks_r_data = (short *) nim_landmarks_r->data;
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
    if (nim_landmarks_r->datatype == NIFTI_TYPE_INT32) {
        int * nim_landmarks_r_data = (int *) nim_landmarks_r->data;
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
    if (nim_layers_r->datatype == NIFTI_TYPE_FLOAT32) {
        float * nim_layers_r_data = (float *) nim_layers_r->data;
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
        short * nim_layers_r_data = (short *) nim_layers_r->data;
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
        int * nim_layers_r_data = (int *) nim_layers_r->data;
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
    cout << "  There are  " << layernumber << " layers  " << endl;

    ////////////////////////////////
    // Allocating necessary files //
    ///////////////////////////////
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

    short * Grow_x_data = (short *) Grow_x->data;
    short * Grow_y_data = (short *) Grow_y->data;
    short * Grow_z_data = (short *) Grow_z->data;

    nifti_image * growfromCenter = nifti_copy_nim_info(nim_layers);
    nifti_image * growfromCenter_thick = nifti_copy_nim_info(nim_layers);

    growfromCenter->datatype = NIFTI_TYPE_INT16;
    growfromCenter_thick->datatype = NIFTI_TYPE_INT16;

    growfromCenter->nbyper = sizeof(short);
    growfromCenter_thick->nbyper = sizeof(short);

    growfromCenter->data = calloc(growfromCenter->nvox, growfromCenter->nbyper);
    growfromCenter_thick->data = calloc(growfromCenter_thick->nvox, growfromCenter_thick->nbyper);

    short * growfromCenter_data = (short *) growfromCenter->data;
    short * growfromCenter_thick_data = (short *) growfromCenter_thick->data;

    ///////////////////////////////
    // Prepare growing variables //
    ///////////////////////////////
    float x1g = 0.;
    float y1g = 0.;
    float z1g = 0.;

    float min_val = 0.;
    float dist_min2 = 0.;
    float dist_i = 0.;
    float dist_p1 = 0.;

    int grow_vinc = 3;
    int grow_vinc_area = 1;
    // int vinc_max = 40 ;
    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);

    ///////////////////////////////////////
    // Growing from Center cross columns //
    ///////////////////////////////////////
    cout << "  Growing from center " << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                if (*(nim_landmarks_data + nxy * iz + nx * ix + iy) == 1 &&  abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber/2)) < 2) {  //defining seed at center landmark
                    *(growfromCenter_data + nxy * iz + nx * ix + iy) = 1.;
                    *(Grow_x_data + nxy * iz + nx * ix + iy) = ix;
                    *(Grow_y_data + nxy * iz + nx * ix + iy) = iy;
                    *(Grow_z_data + nxy * iz + nx * ix + iy) = iz;
                }
            }
        }
    }
    cout << "  Growing " << vinc_max << " iteration " << flush;

    for (int grow_i = 1; grow_i < vinc_max; grow_i++) {
        cout << "\r  " << grow_i << " " << flush;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    if (abs((int) (*(nim_layers_data + nxy * iz + nx * ix + iy) - layernumber/2)) < 2 && *(growfromCenter_data + nxy * iz + nx * ix + iy) == 0 && *(nim_landmarks_data + nxy * iz + nx * ix + iy) < 2) {
                        // Only grow into areas that are GM and that have not
                        // been grown into, yet... and it should stop as soon
                        // as it hits the border
                        for (int iy_i = max(0, iy - grow_vinc_area); iy_i <= min(iy + grow_vinc_area, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - grow_vinc_area); ix_i <= min(ix + grow_vinc_area, sizeRead - 1); ++ix_i) {
                                for (int iz_i = max(0, iz - grow_vinc_area); iz_i <= min(iz + grow_vinc_area, sizeSlice - 1); ++iz_i) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    if (*(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i) == grow_i  && *(nim_landmarks_data + nxy * iz + nx * ix + iy) < 2) {
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
                        if (dist_min2 < 1.7) { // ???? I DONT REMEMBER WHY I NEED THIS ????
                            // distDebug(0, islice, iy, ix) = dist_min2 ;
                            *(growfromCenter_data + nxy * iz + nx * ix + iy) = grow_i+1;
                            *(Grow_x_data + nxy * iz + nx * ix + iy) = *(Grow_x_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_y_data + nxy * iz + nx * ix + iy) = *(Grow_y_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_z_data + nxy * iz + nx * ix + iy) = *(Grow_z_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                        // cout << " ix   " << ix << " iy   " << iy << "    " << * (WMkoord0_data + nxy * islice + nx * (int)x1g + (int)y1g) << endl;
                    }
                }
            }
        }
    }
    cout << endl << "  Growing is done." << flush;

    if (verbose == 1) {
        cout << "  Writing output " << endl;
        const char * fout_5 = "coordinates_1_path.nii";
        if (nifti_set_filenames(growfromCenter, fout_5, 1, 1)) {
            return 1;
        }
        nifti_image_write(growfromCenter);
    }

    /////////////////////
    // Smooth columns  //
    /////////////////////
    // NOTE(Renzo): In the future this should be done only within connected
    // areas. Right now there might be a problem, when the center of two GM.

    // cout << "  Smooth columns " << endl;  // Ribbons is closer than 5 voxels (vinc_sm)

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
    float * smoothed_data = (float *) smoothed->data;
    float * gausweight_data = (float *) gausweight->data;

    // float kernel_size = 10;  // corresponds to one voxel sice.
    int FWHM_val = 1;
    int vinc_sm = 5;// if voxel is too far away, I ignore it.
    dist_i = 0.;
    cout << "    vinc_sm " << vinc_sm << endl;
    cout << "    FWHM_val " << FWHM_val << endl;
    // cout << "    DEBUG " << dist(1., 1., 1., 1., 2., 1., dX, dY, dZ) << endl;
    cout << "    Here 2 " << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                // * (smoothed_data + nxy * iz + nx * ix + iy) = 0 ;
                if (*(growfromCenter_data + nxy * iz + nx * ix + iy)  > 0) {
                    for (int iz_i = max(0, iz - vinc_sm); iz_i <= min(iz + vinc_sm, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc_sm); iy_i <= min(iy + vinc_sm, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc_sm); ix_i <= min(ix + vinc_sm, sizeRead - 1); ++ix_i) {
                                if (*(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i)  > 0) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    // cout << "    Debug  4 " << gaus(dist_i , FWHM_val) << endl;
                                    // cout << "    Debug  5 " << dist_i << endl;
                                    // if (* (nim_input_data + nxy * iz + nx * ix + iy) == 3) cout << "debug  4b " << endl;
                                    // dummy = * (layer_data + nxy * iz_i + nx * ix_i + iy_i);
                                    *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                            }
                        }
                    }
                    if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) {
                        *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy)/ *(gausweight_data + nxy * iz + nx * ix + iy);
                    }
                }
                // if (* (nim_layers_r_data + nxy * iz + nx * ix + iy) <= 0) {
                //     * (smoothed_data + nxy * iz + nx * ix + iy) = * (lateralCoord_data + nxy * iz + nx * ix + iy) ;
                // }
            }
        }
    }
    cout << "    Here 2." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(growfromCenter_data + nxy * iz + nx * ix + iy) > 0) {
                    *(growfromCenter_data + nxy * iz + nx * ix + iy) = (int) *(smoothed_data + nxy * iz + nx * ix + iy);
                }
            }
        }
    }

    if (verbose == 1) {
        const char * fout_6 = "coordinates_2_path_smooth.nii";
        if (nifti_set_filenames(growfromCenter, fout_6, 1, 1)) return 1;
        nifti_image_write(growfromCenter);
    }

    /////////////////////////////////////
    // Extending columns across layers //
    /////////////////////////////////////
    // NOTE(Renzo): This is not perfect yet, because it has only 4 directions
    // to grow thus ther might be orientation biases.
    cout << " extending columns across layers  " << endl;

    nifti_image * hairy_brain = nifti_copy_nim_info(nim_layers);
    hairy_brain->datatype = NIFTI_TYPE_INT16;
    hairy_brain->nbyper = sizeof(short);
    hairy_brain->data = calloc(hairy_brain->nvox, hairy_brain->nbyper);
    short * hairy_brain_data = (short *) hairy_brain->data;
    nifti_image * hairy_brain_dist = nifti_copy_nim_info(nim_layers);
    hairy_brain_dist->datatype = NIFTI_TYPE_FLOAT32;
    hairy_brain_dist->nbyper = sizeof(float);
    hairy_brain_dist->data = calloc(hairy_brain_dist->nvox, hairy_brain_dist->nbyper);
    float * hairy_brain_dist_data = (float *) hairy_brain_dist->data;

    // This is an upper limit of the cortical thickness
    dist_min2 = 10000.;

    // This is the area the algorithm looks for the closest middele layer.
    // The only problem, when this is too big is the longer calculation time.
    int vinc_thickness = 30;

    // This is step size that neigbouring GM voxels need to be to be classified
    // as one side of the GM bank.
    int vinc_steps = 1;

    int cloasest_coord = 0.;
    int there_is_close_noigbour = 0;
    float average_neigbours = 0;
    float average_val = 0;

    cout << "  Growing " << vinc_thickness << " iteration aross layers " << flush;

    int vinc_steps_g = 1;
    int vinc_sm_g = 25;
    int pref_ratio = 0;

    // For estimation of time
    int nvoxels_to_go_across = 0;
    int running_index = 0;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > 1 && *(nim_layers_data + nxy * iz + nx * ix + iy) < layernumber - 1) nvoxels_to_go_across++;
            }
        }
    }

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > 1 && *(nim_layers_data + nxy * iz + nx * ix + iy) < layernumber - 1) {
                    running_index++;
                    if ((running_index * 100) / nvoxels_to_go_across != pref_ratio) {
                        cout << "\r" << (running_index * 100)/nvoxels_to_go_across << "% " << flush;
                        pref_ratio = (running_index * 100)/nvoxels_to_go_across;
                    }
                    /////////////////////////////////////////////////
                    // Find area that is not from the other sulcus //
                    /////////////////////////////////////////////////

                    // Preparation of dummy vicinity dile, resting it with zeros
                    for (int iz_i = max(0, iz - vinc_sm_g - vinc_steps); iz_i <= min(iz + vinc_sm_g + vinc_steps, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc_sm_g - vinc_steps); iy_i <= min(iy + vinc_sm_g + vinc_steps, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc_sm_g - vinc_steps); ix_i <= min(ix + vinc_sm_g + vinc_steps, sizeRead - 1); ++ix_i) {
                                *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = 0;
                            }
                        }
                    }
                    min_val = 0;

                    // Iteration loop that determines a local patch of connected
                    // voxels, exlcudid covels from opposite GM bank.
                    // NOTE(Renzo): This loop takes forever
                    // TODO(Faruk): I need to have a look at this.
                    *(hairy_brain_data + nxy * iz + nx * ix + iy) = 1;
                    for (int K_ = 0; K_ < vinc_sm_g; K_++) {
                        for (int iz_ii = max(0, iz - vinc_sm_g); iz_ii <= min(iz + vinc_sm_g, sizeSlice - 1); ++iz_ii) {
                            for (int iy_ii = max(0, iy - vinc_sm_g); iy_ii <= min(iy + vinc_sm_g, sizePhase - 1); ++iy_ii) {
                                for (int ix_ii = max(0, ix - vinc_sm_g); ix_ii <= min(ix + vinc_sm_g, sizeRead - 1); ++ix_ii) {
                                    if (*(hairy_brain_data + nxy * iz_ii + nx * ix_ii + iy_ii) == 1) {
                                        for (int iz_i = max(0, iz_ii - vinc_steps); iz_i <= min(iz_ii + vinc_steps, sizeSlice - 1); ++iz_i) {
                                            for (int iy_i = max(0, iy_ii - vinc_steps); iy_i <= min(iy_ii + vinc_steps, sizePhase - 1); ++iy_i) {
                                                for (int ix_i = max(0, ix_ii - vinc_steps); ix_i <= min(ix_ii + vinc_steps, sizeRead - 1); ++ix_i) {
                                                    if (dist((float)ix_ii, (float)iy_ii, (float)iz_ii, (float)ix_i, (float)iy_i, (float)iz_i, 1, 1, 1) <= 1 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) > 1 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) < layernumber - 1) {
                                                        *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    // Only grow into areas that are GM and that have not been
                    // grown into, yet... and it should stop as soon as it hits
                    // the border
                    for (int iy_i = max(0, iy - vinc_sm_g); iy_i <= min(iy + vinc_sm_g, sizePhase - 1); ++iy_i) {
                        for (int ix_i = max(0, ix - vinc_sm_g); ix_i <= min(ix + vinc_sm_g, sizeRead - 1); ++ix_i) {
                            for (int iz_i = max(0, iz - vinc_sm_g); iz_i <= min(iz + vinc_sm_g, sizeSlice - 1); ++iz_i) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                if (*(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) == 1 && dist_i < dist_min2 && *(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i) > 0) {
                                    dist_min2 = dist_i;
                                    x1g = ix_i;
                                    y1g = iy_i;
                                    z1g = iz_i;
                                    dist_p1 = dist_min2;
                                    min_val = *(growfromCenter_data + nxy * iz_i + nx * ix_i + iy_i);
                                }
                            }
                        }
                    }
                    *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) = min_val;
                }
            }
        }
    }
    cout << endl;  // to close the online output

    if (verbose == 1) {
        const char * fout_7 = "coordinates_3_ thick.nii";
        if (nifti_set_filenames(growfromCenter_thick, fout_7, 1, 1)) return 1;
        nifti_image_write(growfromCenter_thick);
    }

    ////////////////////////////////////////////////////////////
    // Smooth columns within the thick cortex.                //
    // This smoothing is done to correct for Pytagoras errors //
    // This smoothing is within GM banks only                 //
    // (which makes it slow)                                  //
    ////////////////////////////////////////////////////////////
    cout << "  Smoothing the thick cortex now..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
            }
        }
    }

    FWHM_val = 3;

    // This is small so that neighboring sulci are not affecting each other
    vinc_sm = 8;  // max(1., 2. * FWHM_val/dX);

    cout << "    vinc_sm " << vinc_sm << endl;
    cout << "    FWHM_val " << FWHM_val << endl;
    cout << "    Starting extended now  " << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead-0; ++ix) {
                *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                // * (smoothed_data + nxy * iz + nx * ix + iy) = 0 ;
                if (*(growfromCenter_thick_data + nxy * iz + nx * ix + iy) > 0) {

                    /////////////////////////////////////////////////
                    // Find area that is not from the other sulcus //
                    /////////////////////////////////////////////////

                    // Preparation of dummy vicinity
                    for (int iz_i = max(0, iz - vinc_sm - vinc_steps); iz_i <= min(iz + vinc_sm + vinc_steps, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc_sm - vinc_steps); iy_i <= min(iy + vinc_sm + vinc_steps, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc_sm - vinc_steps); ix_i <= min(ix + vinc_sm + vinc_steps, sizeRead - 1); ++ix_i) {
                                *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = 0;
                            }
                        }
                    }
                    *(hairy_brain_data + nxy * iz + nx * ix + iy) = 1;
                    for (int K_ = 0; K_ < vinc_sm; K_++) {
                        for (int iz_ii = max(0, iz - vinc_sm); iz_ii <= min(iz + vinc_sm, sizeSlice - 1); ++iz_ii) {
                            for (int iy_ii = max(0, iy - vinc_sm); iy_ii <= min(iy + vinc_sm, sizePhase - 1); ++iy_ii) {
                                for (int ix_ii = max(0, ix - vinc_sm); ix_ii <= min(ix + vinc_sm, sizeRead - 1); ++ix_ii) {
                                    if (*(hairy_brain_data + nxy * iz_ii + nx * ix_ii + iy_ii) == 1) {
                                        for (int iz_i = max(0, iz_ii - vinc_steps); iz_i <= min(iz_ii + vinc_steps, sizeSlice - 1); ++iz_i) {
                                            for (int iy_i = max(0, iy_ii - vinc_steps); iy_i <= min(iy_ii + vinc_steps, sizePhase - 1); ++iy_i) {
                                                for (int ix_i = max(0, ix_ii - vinc_steps); ix_i <= min(ix_ii + vinc_steps, sizeRead - 1); ++ix_i) {
                                                    if (dist((float)ix_ii, (float)iy_ii, (float)iz_ii, (float)ix_i, (float)iy_i, (float)iz_i, 1, 1, 1) <= 1 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) > 1 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) < layernumber - 1) {
                                                        *(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Smoothing within each layer and within the local patch
                    int layernumber_i = *(nim_layers_data + nxy * iz + nx * ix + iy);

                    for (int iz_i = max(0, iz - vinc_sm); iz_i <= min(iz + vinc_sm, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc_sm); iy_i <= min(iy + vinc_sm, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc_sm); ix_i <= min(ix + vinc_sm, sizeRead - 1); ++ix_i) {
                                if (*(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) == 1 && abs((int) *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) - layernumber_i) < 2 && *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) > 0) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                    *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                            }
                        }
                    }
                    if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) {
                        *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy)/ *(gausweight_data + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }
    }
    cout << "    Extended now." << endl;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(growfromCenter_thick_data + nxy * iz + nx * ix + iy) > 0) {
                    *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    cout << "    Smoothing done." << endl;

    if (verbose == 1) {
        const char * fout_9 = "coordinates_4_thick_smoothed.nii";
        if (nifti_set_filenames(growfromCenter_thick, fout_9, 1, 1)) {
            return 1;
        }
        nifti_image_write(growfromCenter_thick);
    }

    //////////////////////////
    // Grow final outer rim //
    //////////////////////////

    // Note(Renzo): I could not fill this earlier because it would have resulted
    // in leakage from the opposite back. Thus I always left a safety corridor.
    // Which is filled in now.

    int vinc_rim = 2;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                *(hairy_brain_data + nxy * iz + nx * ix + iy) = *(growfromCenter_thick_data + nxy * iz + nx * ix + iy);
            }
        }
    }
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > 0 && *(growfromCenter_thick_data + nxy * iz + nx * ix + iy) == 0) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    min_val = 0;
                    // Only grow into areas that are GM and that have not been
                    // grown into, yet... And it should stop as soon as it hits
                    // the border.
                    for (int iy_i = max(0, iy - vinc_rim); iy_i <= min(iy + vinc_rim, sizePhase - 1); ++iy_i) {
                        for (int ix_i = max(0, ix - vinc_rim); ix_i <= min(ix + vinc_rim, sizeRead - 1); ++ix_i) {
                            for (int iz_i = max(0, iz - vinc_rim); iz_i <= min(iz + vinc_rim, sizeSlice - 1); ++iz_i) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                if (dist_i < dist_min2 && *(nim_layers_data + nxy * iz_i + nx * ix_i + iy_i) > 0 && *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i) > 0) {
                                    dist_min2 = dist_i;
                                    x1g = ix_i;
                                    y1g = iy_i;
                                    z1g = iz_i;
                                    dist_p1 = dist_min2;
                                    min_val = *(growfromCenter_thick_data + nxy * iz_i + nx * ix_i + iy_i);
                                }
                            }
                        }
                    }
                    *(hairy_brain_data + nxy * iz + nx * ix + iy) = min_val;
                }
            }
        }
    }

    if (verbose == 1) {
        const char * fout_3 = "coordinates_5_extended.nii";
        if (nifti_set_filenames(hairy_brain, fout_3, 1, 1)) return 1;
        nifti_image_write(hairy_brain);
    }

    ////////////////////////////////////////////////////////////
    // Resample the number of columns to the user given value //
    ////////////////////////////////////////////////////////////
    if (Ncolumns >0)  {
        cout << "   Resampling the number of columns..." << endl;
        int max_columns = 0;
        int min_columns = 100000000;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    if (*(hairy_brain_data + nxy * iz + nx * ix + iy) >  0) {
                        if ((int) *(hairy_brain_data + nxy * iz + nx * ix + iy) >  max_columns) max_columns = (int) *(hairy_brain_data + nxy * iz + nx * ix + iy);
                        if ((int) *(hairy_brain_data + nxy * iz + nx * ix + iy) < min_columns) min_columns = (int) *(hairy_brain_data + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }

        cout << "   Max = " << max_columns << " | Min = " << min_columns << endl;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    if (*(hairy_brain_data + nxy * iz + nx * ix + iy) >  0) {
                        *(hairy_brain_data + nxy * iz + nx * ix + iy) = (*(hairy_brain_data + nxy * iz + nx * ix + iy) - (short)min_columns) * (short)(Ncolumns - 1)/(short)(max_columns-min_columns) +1;
                    }
                }
            }
        }
    }

    const char * fout_4 = "coordinates_final.nii";
    if (nifti_set_filenames(hairy_brain, fout_4, 1, 1)) {
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
                + (z1 - z2) * (z1-z2) * dZ * dZ);
}


float gaus(float distance, float sigma) {
    return (1./(sigma * sqrt(2. * 3.141592))
            * exp(-0.5 * distance * distance / (sigma * sigma)));
}
