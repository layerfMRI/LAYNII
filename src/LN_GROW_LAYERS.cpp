

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_GROW_LAYERS: Example of gray matter layering.\n"
    "\n"
    "    This program calculates the layers based on GM and CSF border line \n"
    "    rim files.\n"
    "\n"
    "Usage:\n"
    "    LN_GROW_LAYERS -rim rim.nii -N 21 \n"
    "    LN_GROW_LAYERS -rim rim.nii -N 21 -vinc 40 \n"
    " \n"
    "Test usage in the test_data folder: ../LN_GROW_LAYERS -rim sc_rim.nii \n"
    "\n"
    "Example application in a blog post: https://layerfmri.com/quick-layering/ \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -rim        : Specify input dataset.\n"
    "    -vinc       : Size of vicinity. Default is 40. This is the maximum\n"
    "                  thickness of the cortex in units of voxels. Smaller\n"
    "                  the number the faster the program.\n"
    "    -N          : (Optional) Number of layers. Default is 20.\n"
    "                  In visual cortex you might want to use less. \n"
    "                  Maximum accuracy is 1/100 for now.\n"
    "    -thin       : (Optional) When the distance between layers is less\n"
    "                  than the voxel thickness. This deals with missing\n"
    "                  layers next to the inner most and outer most layers.\n"
    "    -threeD     : Do layer calculations in 3D. Default is 2D.\n"
    "    -debug      : If you want to see the growing of the respective\n"
    "                  tissue types, it is writen out.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image* nim_input_i = NULL;
    char* fin = NULL;
    int ac, Nlayer_real = 20, vinc_int = 50;
    int threeD = 0, thinn_option = 0, centroid_option = 0, debug = 0;
    if (argc< 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac< argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-N")) {
            if (++ac >= argc) {
                return 1;
            }
            Nlayer_real = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-vinc")) {
            if (++ac >= argc) {
                return 1;
            }
            vinc_int = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-threeD")) {
            fprintf(stderr, "Layer calculation will be done in 3D.\n ");
            threeD = 1;
        } else if (!strcmp(argv[ac], "-debug")) {
            fprintf(stderr, "Writing out growing fields.\n ");
            debug = 1;
        } else if (!strcmp(argv[ac], "-thin")) {
            fprintf(stderr, "Correct for extra thin layers.\n ");
            thinn_option = 1;
        } else if (!strcmp(argv[ac], "-centroid")) {
            fprintf(stderr, "Write out another file with centroid depth.\n ");
            centroid_option = 1;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) { fprintf(stderr, "** missing option '-rim'\n");  return 1; }
    // read input dataset, including data
    nim_input_i = nifti_image_read(fin, 1);
    if (!nim_input_i) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_GROW_LAYERS");
    log_nifti_descriptives(nim_input_i);

    // Get dimensions of input
    int size_z = nim_input_i->nz;
    int size_x = nim_input_i->nx;
    int size_y = nim_input_i->ny;
    int size_time = nim_input_i->nt;
    int nx = nim_input_i->nx;
    int nxy = nim_input_i->nx * nim_input_i->ny;
    int nxyz = nim_input_i->nx * nim_input_i->ny * nim_input_i->nz;
    float dX = nim_input_i->pixdim[1];
    float dY = nim_input_i->pixdim[2];
    float dZ = nim_input_i->pixdim[3];

    cout << "  Using " << Nlayer_real << " layers." << endl;
    if (Nlayer_real > 1000) {
        cout << endl << "  Stop! Too many layers (>1000). I can't let you do that." << endl;
    }
    int Nlayer = 1000;  // This is an interim number that will be scaled later

    if (vinc_int != 50) {
        cout << "  Calculate layers up to cortical thicknesses of " << vinc_int << " voxels " << endl;
    }

    // ========================================================================
    // Fix data type issues
    nifti_image* nim_input = copy_nifti_as_int32(nim_input_i);
    int32_t* nim_input_data = static_cast<int32_t*>(nim_input->data);

    // Allocate new nifti images
    nifti_image * nii_layers = copy_nifti_as_int32(nim_input);
    int32_t* nii_layers_data = static_cast<int32_t*>(nii_layers->data);

    // ========================================================================
    ///////////////////////////////////////////////////
    // I am doing a the layer calculation in 2D here //
    ///////////////////////////////////////////////////
    if (threeD == 0) {
        nifti_image* growfromWM0 = copy_nifti_as_float32(nim_input);
        float* growfromWM0_data = static_cast<float*>(growfromWM0->data);
        nifti_image * growfromWM1 = copy_nifti_as_float32(nim_input);
        float* growfromWM1_data = static_cast<float*>(growfromWM1->data);

        nifti_image* WMkoord0 = copy_nifti_as_int16(nim_input);
        nifti_image* WMkoord1 = copy_nifti_as_int16(nim_input);
        nifti_image* WMkoord2 = copy_nifti_as_int16(nim_input);
        nifti_image* WMkoord3 = copy_nifti_as_int16(nim_input);

        int16_t* WMkoord0_data = static_cast<int16_t*>(WMkoord0->data);
        int16_t* WMkoord1_data = static_cast<int16_t*>(WMkoord1->data);
        int16_t* WMkoord2_data = static_cast<int16_t*>(WMkoord2->data);
        int16_t* WMkoord3_data = static_cast<int16_t*>(WMkoord3->data);


        nifti_image* growfromGM0 = copy_nifti_as_float32(nim_input);
        nifti_image* growfromGM1 = copy_nifti_as_float32(nim_input);
        float* growfromGM0_data = static_cast<float*>(growfromGM0->data);
        float* growfromGM1_data = static_cast<float*>(growfromGM1->data);

        nifti_image* GMkoord0 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoord1 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoord2 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoord3 = copy_nifti_as_int16(nim_input);

        int16_t* GMkoord0_data = static_cast<int16_t*>(GMkoord0->data);
        int16_t* GMkoord1_data = static_cast<int16_t*>(GMkoord1->data);
        int16_t* GMkoord2_data = static_cast<int16_t*>(GMkoord2->data);
        int16_t* GMkoord3_data = static_cast<int16_t*>(GMkoord3->data);
        // --------------------------------------------------------------------
        // Coordinates
        float x1g = 0., y1g = 0.;
        float x2g = 0., y2g = 0.;
        float x3g = 0., y3g = 0.;

        cout << "  Until here 2 " << endl;
        // Reduce mask to contain only areas close to the surface.
        cout << "  Select GM regions..." << endl;

        // This is the distance from every voxel the algorithm is applied on
        // Just to make it faster and not loop over all voxels.
        int vinc = vinc_int;

        float dist_i = 0.;
        float dist_min = 0.;
        float dist_min1 = 0.;
        float dist_min2 = 0.;
        float dist_min3 = 0.;
        float dist_max = 0.;
        float dist_p1 = 0.;

        cout << "  Start growing from WM..." << endl;

        // Setting zero
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromWM0_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromWM1_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromGM0_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromGM1_data + nxy * iz + nx * ix + iy) = 0;
                }
            }
        }

        //////////////////
        // Grow from WM //
        //////////////////
        int grow_vinc = 2;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2) {
                        *(growfromWM0_data + nxy * iz + nx * ix + iy) = 1.;
                        *(WMkoord0_data + nxy * iz + nx * ix + iy) = ix;
                        *(WMkoord1_data + nxy * iz + nx * ix + iy) = iy;
                    }
                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3  && *(growfromWM0_data + nxy * iz + nx * ix + iy) == 0) {
                            // cout << " true " << * (growfromWM0_data + nxy * iz + nx * ix + iy) << endl;
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, size_y); ++ix_i) {
                                    if (*(growfromWM0_data + nxy * iz + nx * ix_i + iy_i) == (float)grow_i) {
                                        dist_i = dist2d((float)ix, (float)iy, (float)ix_i, (float)iy_i);
                                        if (dist_i< dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                            if (dist_min2< 1.4) {
                                // distDebug(0, iz, iy, ix) = dist_min2 ;
                                *(growfromWM0_data + nxy * iz + nx * ix + iy) = (float)grow_i+1;
                                *(WMkoord0_data + nxy * iz + nx * ix + iy) = *(WMkoord0_data + nxy * iz + nx * (int)x1g + (int)y1g);
                                *(WMkoord1_data + nxy * iz + nx * ix + iy) = *(WMkoord1_data + nxy * iz + nx * (int)x1g + (int)y1g);
                            }
                            // cout << " ix " << ix << " iy " << iy << " " << * (WMkoord0_data + nxy * iz + nx * (int)x1g + (int)y1g) << endl;
                        }
                    }
                }
            }
        }

        ///////////////////
        // grow from CSF //
        ///////////////////
        cout << "  Start growing from CSF..." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                        *(growfromGM0_data + nxy * iz + nx * ix + iy) = 1.;
                        *(GMkoord0_data + nxy * iz + nx * ix + iy) = ix;
                        *(GMkoord1_data + nxy * iz + nx * ix + iy) = iy;
                    }
                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        dist_min2 = 10000.;
                        x1g = 0.;
                        y1g = 0;
                        if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3  && *(growfromGM0_data + nxy * iz + nx * ix + iy) == 0) {
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, size_y); ++ix_i) {
                                    if (*(growfromGM0_data + nxy * iz + nx * ix_i + iy_i) == (float)grow_i) {
                                        dist_i = dist2d((float)ix, (float)iy, (float)ix_i, (float)iy_i);
                                        if (dist_i< dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                            if (dist_min2< 1.4) {
                                *(growfromGM0_data + nxy * iz + nx * ix + iy) = (float)grow_i+1;

                                *(GMkoord0_data + nxy * iz + nx * ix + iy) = *(GMkoord0_data + nxy * iz + nx * (int)x1g + (int)y1g);
                                *(GMkoord1_data + nxy * iz + nx * ix + iy) = *(GMkoord1_data + nxy * iz + nx * (int)x1g + (int)y1g);
                            }
                        }
                    }
                }
            }
        }

        if (debug > 0) {
            save_output_nifti(fin, "debug_WM_pre_pythagoras", growfromWM1, false);
            save_output_nifti(fin, "debug_GM_pre_pythagoras", growfromGM1, false);
        }

        ////////////////////////////////////////////////////
        // Wabble across neighbouring voxels of closest WM //
        // to account for Pytagoras errors /////////////////
        ////////////////////////////////////////////////////
        cout << "  Correcting for Pytagoras error..." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromWM1_data + nxy * iz + nx * ix + iy) = *(growfromWM0_data + nxy * iz + nx * ix + iy);
                    *(WMkoord2_data + nxy * iz + nx * ix + iy) = *(WMkoord0_data + nxy * iz + nx * ix + iy);
                    *(WMkoord3_data + nxy * iz + nx * ix + iy) = *(WMkoord1_data + nxy * iz + nx * ix + iy);

                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        if (*(WMkoord1_data + nxy * iz + nx * ix + iy) != 0) {
                            dist_min2 = 10000.;
                            x1g = 0;
                            y1g = 0;

                            for (int iy_i = max(0, (*(WMkoord3_data + nxy * iz + nx * ix + iy)) - grow_vinc); iy_i < min((*(WMkoord3_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, (*(WMkoord2_data + nxy * iz + nx * ix + iy)) - grow_vinc); ix_i < min((*(WMkoord2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_y); ++ix_i) {
                                    if (*(nim_input_data + nxy * iz + nx * ix_i + iy_i) == 2) {

                                        dist_i = dist2d((float)ix, (float)iy, (float)ix_i, (float)iy_i);
                                        if (dist_i< dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                            *(growfromWM1_data + nxy * iz + nx * ix + iy) = dist2d((float)ix, (float)iy, (float)x1g, (float)y1g);
                            *(WMkoord2_data + nxy * iz + nx * ix + iy) = *(WMkoord2_data + nxy * iz + nx * (int)x1g + (int)y1g);
                            *(WMkoord3_data + nxy * iz + nx * ix + iy) = *(WMkoord3_data + nxy * iz + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
        // cout << "  Running until here..." << endl;

        ////////////////////////////////////////////////////
        // Wabble across neighbouring voxels of closest GM //
        // to account for Pytagoras errors /////////////////
        ////////////////////////////////////////////////////
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromGM1_data + nxy * iz + nx * ix + iy) = *(growfromGM0_data + nxy * iz + nx * ix + iy);
                    *(GMkoord2_data + nxy * iz + nx * ix + iy) = *(GMkoord0_data + nxy * iz + nx * ix + iy);
                    *(GMkoord3_data + nxy * iz + nx * ix + iy) = *(GMkoord1_data + nxy * iz + nx * ix + iy);
                }
            }
        }
        // cout << "  Running also until here..." << endl;
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        if (*(GMkoord1_data + nxy * iz + nx * ix + iy) != 0) {
                            dist_min2 = 10000.;
                            x1g = 0;
                            y1g = 0;
                            for (int iy_i = max(0, (*(GMkoord3_data + nxy * iz + nx * ix + iy)) - grow_vinc); iy_i < min((*(GMkoord3_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, (*(GMkoord2_data + nxy * iz + nx * ix + iy)) - grow_vinc); ix_i < min((*(GMkoord2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_y); ++ix_i) {
                                    if (*(nim_input_data + nxy * iz + nx * ix_i + iy_i) == 1) {
                                        dist_i = dist2d((float)ix, (float)iy, (float)ix_i, (float)iy_i);
                                        if (dist_i< dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = ix_i;
                                            y1g = iy_i;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                            *(growfromGM1_data + nxy * iz + nx * ix + iy) = dist2d((float)ix, (float)iy, (float)x1g, (float)y1g);
                            *(GMkoord2_data + nxy * iz + nx * ix + iy) = *(GMkoord2_data + nxy * iz + nx * (int)x1g + (int)y1g);
                            *(GMkoord3_data + nxy * iz + nx * ix + iy) = *(GMkoord3_data + nxy * iz + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
        // cout << "  Running also until here 3..." << endl;

        int GMK2_i, GMK3_i, WMK2_i, WMK3_i;
        float GMK2_f, GMK3_f, WMK2_f, WMK3_f, ix_f, iy_f;

        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) {
                        // Equi_dist_layers(0, iz, iy, ix) = 19 * (1- dist((float)ix, (float)iy, (float)GMkoord(2, iz, iy, ix) , (float)GMkoord(3, iz, iy, ix) )/ (dist((float)ix, (float)iy, (float)GMkoord(2, iz, iy, ix)    , (float)GMkoord(3, iz, iy, ix) ) + dist((float)ix, (float)iy, (float)WMkoord(2, iz, iy, ix)        , (float)WMkoord(3, iz, iy, ix)))) + 2  ;

                        GMK2_i = *(GMkoord2_data + nxy * iz + nx * ix + iy);
                        GMK3_i = *(GMkoord3_data + nxy * iz + nx * ix + iy);
                        WMK2_i = *(WMkoord2_data + nxy * iz + nx * ix + iy);
                        WMK3_i = *(WMkoord3_data + nxy * iz + nx * ix + iy);
                        GMK2_f = (float)GMK2_i;
                        GMK3_f = (float)GMK3_i;
                        WMK2_f = (float)WMK2_i;
                        WMK3_f = (float)WMK3_i;
                        ix_f = (float)ix;
                        iy_f = (float)iy;

                        // cout << " rix_f, iy_f, GMK2_f, GMK3_f " << " " << ix_f << " " << iy_f << " " << GMK2_f << " " << * (GMkoord2_data + nxy * iz + nx * ix + iy) << endl;
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = (Nlayer-1) * (1- dist2d((float)ix, (float)iy, GMK2_f, GMK3_f)/ (dist2d((float)ix, (float)iy, GMK2_f, GMK3_f) + dist2d((float)ix, (float)iy, WMK2_f, WMK3_f))) + 2;
                        // * (nii_layers_data + nxy * iz + nx * ix + iy) = 100 * dist(ix_f, iy_f, GMK2_f, GMK3_f) ;
                    }
                }
            }
        }

        if (debug > 0) {
            save_output_nifti(fin, "debug_WM", growfromWM1, false);
            save_output_nifti(fin, "debug_GM", growfromGM1, false);
        }
    }  // 2D layer calculation is closed.

    //////////////////////////////////////////////
    // Doing a the layer calculation in 3D here //
    //////////////////////////////////////////////
    if (threeD == 1) {
        cout << "  Starting 3D loop..." << endl;
        nifti_image* growfromWM0 = copy_nifti_as_float32(nim_input);
        nifti_image* growfromWM1 = copy_nifti_as_float32(nim_input);
        float* growfromWM0_data = static_cast<float*>(growfromWM0->data);
        float* growfromWM1_data = static_cast<float*>(growfromWM1->data);

        nifti_image * WMkoordx1 = copy_nifti_as_int16(nim_input);
        nifti_image * WMkoordy1 = copy_nifti_as_int16(nim_input);
        nifti_image * WMkoordz1 = copy_nifti_as_int16(nim_input);
        nifti_image * WMkoordx2 = copy_nifti_as_int16(nim_input);
        nifti_image * WMkoordy2 = copy_nifti_as_int16(nim_input);
        nifti_image * WMkoordz2 = copy_nifti_as_int16(nim_input);

        int16_t* WMkoordx1_data = static_cast<int16_t*>(WMkoordx1->data);
        int16_t* WMkoordy1_data = static_cast<int16_t*>(WMkoordy1->data);
        int16_t* WMkoordz1_data = static_cast<int16_t*>(WMkoordz1->data);
        int16_t* WMkoordx2_data = static_cast<int16_t*>(WMkoordx2->data);
        int16_t* WMkoordy2_data = static_cast<int16_t*>(WMkoordy2->data);
        int16_t* WMkoordz2_data = static_cast<int16_t*>(WMkoordz2->data);

        nifti_image* growfromGM0 = copy_nifti_as_float32(nim_input);
        nifti_image* growfromGM1 = copy_nifti_as_float32(nim_input);
        float* growfromGM0_data = static_cast<float*>(growfromGM0->data);
        float* growfromGM1_data = static_cast<float*>(growfromGM1->data);

        nifti_image* GMkoordx1 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoordy1 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoordz1 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoordx2 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoordy2 = copy_nifti_as_int16(nim_input);
        nifti_image* GMkoordz2 = copy_nifti_as_int16(nim_input);

        int16_t* GMkoordx1_data = static_cast<int16_t*>(GMkoordx1->data);
        int16_t* GMkoordy1_data = static_cast<int16_t*>(GMkoordy1->data);
        int16_t* GMkoordz1_data = static_cast<int16_t*>(GMkoordz1->data);
        int16_t* GMkoordx2_data = static_cast<int16_t*>(GMkoordx2->data);
        int16_t* GMkoordy2_data = static_cast<int16_t*>(GMkoordy2->data);
        int16_t* GMkoordz2_data = static_cast<int16_t*>(GMkoordz2->data);

        // Coordinates
        float x1g = 0., y1g = 0., z1g = 0.;
        float x2g = 0., y2g = 0., z2g = 0.;
        float x3g = 0., y3g = 0., z3g = 0.;

        // Reduce mask to contain only areas close to the surface.
        cout << "  Select GM regions..." << endl;

        // This is the distance from every voxel the algorithm is applied on.
        // Just to make it faster and not loop over all voxels.
        int vinc = vinc_int;

        float dist_i = 0.;
        float dist_min = 0.;
        float dist_min1 = 0.;
        float dist_min2 = 0.;
        float dist_min3 = 0.;
        float dist_max = 0.;
        float dist_p1 = 0.;

        cout << "  Start growing from WM..." << endl;

        // Setting zero
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromWM0_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromWM1_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromGM0_data + nxy * iz + nx * ix + iy) = 0;
                    *(growfromGM1_data + nxy * iz + nx * ix + iy) = 0;
                }
            }
        }

        /////////////////////////
        // Closing 3D surfaces //
        /////////////////////////
        // cout << "  Closing surfaces..." << endl;
        // int isWMb, isCSFb, isb;
        // for (int iterations = 0; iterations < 100; ++iterations) {
        //     for (int iz = 0; iz < size_z; ++iz) {
        //         for (int iy = 0; iy < size_x; ++iy) {
        //             for (int ix = 0; ix < size_y; ++ix) {
        //                 if (*(nim_input_data + nxy * iz + nx * ix + iy) == 0) {
        //                     isWMb = 0;
        //                     isCSFb = 0;
        //                     isb = 0;
        //                     for (int iz_i = max(0, iz-2); iz_i < min(size_z, iz+2); ++iz_i) {
        //                         for (int iy_i = max(0, iy-2); iy_i < min(size_x, iy+2); ++iy_i) {
        //                             for (int ix_i = max(0, ix-2); ix_i < min(size_y, ix+2); ++ix_i) {
        //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 3) isb = 1;
        //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 1) isCSFb = 1;
        //                                 if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 2) isWMb = 1;
        //                             }
        //                         }
        //                     }
        //                     if (isb == 1 && isCSFb == 1) *(nim_input_data + nxy * iz + nx * ix + iy) = 1;
        //                     if (isb == 1 &&  isWMb == 1) *(nim_input_data + nxy * iz + nx * ix + iy) = 2;
        //                     if (isCSFb == 1 && isWMb == 1) *(nim_input_data + nxy * iz + nx * ix + iy) = 0;// cout << " THIS BRAIN IS WIERD " << endl ;
        //
        //
        //                 }
        //             }
        //         }
        //     }
        // }

        //////////////////
        // Grow from WM //
        //////////////////
        int grow_vinc = 2;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2) {
                        *(growfromWM0_data + nxy * iz + nx * ix + iy) = 1.;
                        *(WMkoordx1_data + nxy * iz + nx * ix + iy) = ix;
                        *(WMkoordy1_data + nxy * iz + nx * ix + iy) = iy;
                        *(WMkoordz1_data + nxy * iz + nx * ix + iy) = iz;
                    }
                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        dist_min2 = 10000.;
                        x1g = 0;
                        y1g = 0;
                        z1g = 0;
                        if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3  && *(growfromWM0_data + nxy * iz + nx * ix + iy) == 0) {
                            // cout << " true " << * (growfromWM0_data + nxy * iz + nx * ix + iy) << endl;
                            for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, size_y); ++ix_i) {
                                    for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc + 1, size_y); ++iz_i) {
                                        if (*(growfromWM0_data + nxy * iz_i + nx * ix_i + iy_i) == (float)grow_i) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            if (dist_i< dist_min2) {
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
                            if (dist_min2< 1.4) {
                                // distDebug(0, iz, iy, ix) = dist_min2 ;
                                *(growfromWM0_data + nxy * iz + nx * ix + iy) = (float)grow_i+1;
                                *(WMkoordx1_data + nxy * iz + nx * ix + iy) = *(WMkoordx1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                                *(WMkoordy1_data + nxy * iz + nx * ix + iy) = *(WMkoordy1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                                *(WMkoordz1_data + nxy * iz + nx * ix + iy) = *(WMkoordz1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            }
                            // cout << " ix " << ix << " iy " << iy << " " << * (WMkoordx1_data + nxy * iz + nx * (int)x1g + (int)y1g) << endl;
                        }
                    }
                }
            }
        }

        ///////////////////
        // Grow from CSF //
        ///////////////////
        cout << "  Start growing from CSF..." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                        *(growfromGM0_data + nxy * iz + nx * ix + iy) = 1.;
                        *(GMkoordx1_data + nxy * iz + nx * ix + iy) = ix;
                        *(GMkoordy1_data + nxy * iz + nx * ix + iy) = iy;
                        *(GMkoordz1_data + nxy * iz + nx * ix + iy) = iz;
                    }
                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        dist_min2 = 10000.;
                        x1g = 0.;
                        y1g = 0;
                        z1g = 0;
                        if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3  && *(growfromGM0_data + nxy * iz + nx * ix + iy) == 0) {
                            for (int iz_i = max(0, iz - grow_vinc); iz_i < min(iz + grow_vinc + 1, size_y); ++iz_i) {
                                for (int iy_i = max(0, iy - grow_vinc); iy_i < min(iy + grow_vinc + 1, size_x); ++iy_i) {
                                    for (int ix_i = max(0, ix - grow_vinc); ix_i < min(ix + grow_vinc + 1, size_y); ++ix_i) {
                                        if (*(growfromGM0_data + nxy * iz_i + nx * ix_i + iy_i) == (float)grow_i) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            if (dist_i< dist_min2) {
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
                                *(growfromGM0_data + nxy * iz + nx * ix + iy) = (float)grow_i+1;
                                *(GMkoordx1_data + nxy * iz + nx * ix + iy) = *(GMkoordx1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                                *(GMkoordy1_data + nxy * iz + nx * ix + iy) = *(GMkoordy1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                                *(GMkoordz1_data + nxy * iz + nx * ix + iy) = *(GMkoordz1_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            }
                        }
                    }
                }
            }
        }

        /////////////////////////////////////////////////////
        // Wabble accross neighbouring voxels of closest WM //
        // to account for Pytagoras errors //////////////////
        /////////////////////////////////////////////////////
        cout << "  Correcting for Pytagoras error..." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromWM1_data + nxy * iz + nx * ix + iy) = *(growfromWM0_data + nxy * iz + nx * ix + iy);
                    *(WMkoordx2_data + nxy * iz + nx * ix + iy) = *(WMkoordx1_data + nxy * iz + nx * ix + iy);
                    *(WMkoordy2_data + nxy * iz + nx * ix + iy) = *(WMkoordy1_data + nxy * iz + nx * ix + iy);
                    *(WMkoordz2_data + nxy * iz + nx * ix + iy) = *(WMkoordz1_data + nxy * iz + nx * ix + iy);
                }
            }
        }
        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        if (*(WMkoordy1_data + nxy * iz + nx * ix + iy) != 0) {
                            dist_min2 = 10000.;
                            x1g = 0;
                            y1g = 0;
                            z1g = 0;
                            for (int iy_i = max(0, (*(WMkoordy2_data + nxy * iz + nx * ix + iy)) - grow_vinc); iy_i < min((*(WMkoordy2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, (*(WMkoordx2_data + nxy * iz + nx * ix + iy)) - grow_vinc); ix_i < min((*(WMkoordx2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_y); ++ix_i) {
                                    for (int iz_i = max(0, (*(WMkoordz2_data + nxy * iz + nx * ix + iy)) - grow_vinc); iz_i < min((*(WMkoordz2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_z); ++iz_i) {
                                        if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 2) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            if (dist_i< dist_min2) {
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
                            *(growfromWM1_data + nxy * iz + nx * ix + iy) = dist((float)ix, (float)iy, (float)iz, (float)x1g, (float)y1g, (float)z1g, dX, dY, dZ);
                            *(WMkoordx2_data + nxy * iz + nx * ix + iy) = *(WMkoordx2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(WMkoordy2_data + nxy * iz + nx * ix + iy) = *(WMkoordy2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(WMkoordz2_data + nxy * iz + nx * ix + iy) = *(WMkoordz2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
        cout << "  Running until here..." << endl;

        ///////////////////////////////////////////////////////
        // Wabble accross neighbouring voxeles of closest GM //
        // To account for Pytagoras errors ////////////////////
        ///////////////////////////////////////////////////////
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(growfromGM1_data + nxy * iz + nx * ix + iy) = *(growfromGM0_data + nxy * iz + nx * ix + iy);
                    *(GMkoordx2_data + nxy * iz + nx * ix + iy) = *(GMkoordx1_data + nxy * iz + nx * ix + iy);
                    *(GMkoordy2_data + nxy * iz + nx * ix + iy) = *(GMkoordy1_data + nxy * iz + nx * ix + iy);
                    *(GMkoordz2_data + nxy * iz + nx * ix + iy) = *(GMkoordz1_data + nxy * iz + nx * ix + iy);
                }
            }
        }
        cout << "  Running also until here..." << endl;

        for (int grow_i = 1; grow_i< vinc; grow_i++) {
            for (int iz = 0; iz < size_z; ++iz) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        if (*(GMkoordy1_data + nxy * iz + nx * ix + iy) != 0) {
                            dist_min2 = 10000.;
                            x1g = 0;
                            y1g = 0;
                            z1g = 0;
                            for (int iy_i = max(0, (*(GMkoordy2_data + nxy * iz + nx * ix + iy)) - grow_vinc); iy_i < min((*(GMkoordy2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_x); ++iy_i) {
                                for (int ix_i = max(0, (*(GMkoordx2_data + nxy * iz + nx * ix + iy)) - grow_vinc); ix_i < min((*(GMkoordx2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_y); ++ix_i) {
                                    for (int iz_i = max(0, (*(GMkoordz2_data + nxy * iz + nx * ix + iy)) - grow_vinc); iz_i < min((*(GMkoordz2_data + nxy * iz + nx * ix + iy)) + grow_vinc + 1, size_y); ++iz_i) {
                                        if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i) == 1) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            if (dist_i< dist_min2) {
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
                            *(growfromGM1_data+ nxy * iz + nx * ix + iy) = dist((float)ix, (float)iy, (float)iz, (float)x1g, (float)y1g, (float)z1g, dX, dY, dZ);
                            *(GMkoordx2_data + nxy * iz + nx * ix + iy) = *(GMkoordx2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(GMkoordy2_data + nxy * iz + nx * ix + iy) = *(GMkoordy2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(GMkoordz2_data + nxy * iz + nx * ix + iy) = *(GMkoordz2_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
        cout << "  Running also until here 3..." << endl;

        int GMK2_i, GMKz2_i, GMK3_i, WMK2_i, WMKz2_i, WMK3_i;
        float GMK2_f, GMKz2_f, GMK3_f, WMK2_f, WMKz2_f, WMK3_f, ix_f, iy_f, iz_f;

        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) {
                        // Equi_dist_layers(0, iz, iy, ix) = 19 * (1- dist((float)ix, (float)iy, (float)GMkoord(2, iz, iy, ix) , (float)GMkoord(3, iz, iy, ix) )/ (dist((float)ix, (float)iy, (float)GMkoord(2, iz, iy, ix)    , (float)GMkoord(3, iz, iy, ix) ) + dist((float)ix, (float)iy, (float)WMkoord(2, iz, iy, ix)        , (float)WMkoord(3, iz, iy, ix)))) + 2  ;

                        GMK2_i = *(GMkoordx2_data + nxy * iz + nx * ix + iy);
                        GMK3_i = *(GMkoordy2_data + nxy * iz + nx * ix + iy);
                        GMKz2_i = *(GMkoordz2_data+ nxy * iz + nx * ix + iy);

                        WMK2_i = *(WMkoordx2_data + nxy * iz + nx * ix + iy);
                        WMK3_i = *(WMkoordy2_data + nxy * iz + nx * ix + iy);
                        WMKz2_i = *(WMkoordz2_data+ nxy * iz + nx * ix + iy);

                        GMK2_f = (float)GMK2_i;
                        GMK3_f = (float)GMK3_i;
                        GMKz2_f = (float)GMKz2_i;

                        WMK2_f = (float)WMK2_i;
                        WMK3_f = (float)WMK3_i;
                        WMKz2_f = (float)WMKz2_i;

                        ix_f = (float)ix;
                        iy_f = (float)iy;
                        iz_f = (float)iz;

                        // cout << " rix_f, iy_f, GMK2_f, GMK3_f " << " " << ix_f << " " << iy_f << " " << GMK2_f << " " << * (GMkoordx2_data + nxy * iz + nx * ix + iy) << endl;
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = (Nlayer-1) * (1- dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ)/ (dist((float)ix, (float)iy, (float)iz, GMK2_f, GMK3_f, GMKz2_f, dX, dY, dZ) + dist((float)ix, (float)iy, (float)iz, WMK2_f, WMK3_f, WMKz2_f, dX, dY, dZ))) + 2;
                        // * (nii_layers_data + nxy * iz + nx * ix + iy) = 100 * dist(ix_f, iy_f, GMK2_f, GMK3_f) ;
                    }
                }
            }
        }
        if (debug > 0) {
            save_output_nifti(fin, "debug_WM", growfromWM1, false);
            save_output_nifti(fin, "debug_GM", growfromWM1, false);
        }
    }  // 3D layer calculation is closed.

    ///////////////////////////////////////////////////////////////
    //// Cleaning negative layers and layers of more than cutoff //
    ///////////////////////////////////////////////////////////////
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_x; ++iy) {
            for (int ix = 0; ix < size_y; ++ix) {
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1 && *(nii_layers_data + nxy * iz + nx * ix + iy) == 0) {
                    *(nii_layers_data + nxy * iz + nx * ix + iy) = Nlayer+1;
                }
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2 && *(nii_layers_data + nxy * iz + nx * ix + iy) == 0) {
                    *(nii_layers_data + nxy * iz + nx * ix + iy) = 1;
                }
            }
        }
    }

    ///////////////////////////////////////////
    // Write out centroid location of layers //
    ///////////////////////////////////////////
    if (centroid_option == 1) {
        nifti_image* centroid = copy_nifti_as_float32(nim_input);
        float* centroid_data = static_cast<float*>(centroid->data);
        float coord = 0.0;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) != 0) {
                        coord = (float)((double) *(nii_layers_data + nxy * iz + nx * ix + iy) /((double)(Nlayer)));
                        if (coord >= 0 && coord<= 1) {  // Taking care of lfoar rounding errors
                            *(centroid_data + nxy * iz + nx * ix + iy) = coord;
                        }
                        if (coord< 0) {  // Taking care of lfoar rounding errors
                            *(centroid_data + nxy * iz + nx * ix + iy) = 0.0;
                        }
                        if (coord > 1) {  // Taking care of lfoar rounding errors
                            *(centroid_data + nxy * iz + nx * ix + iy) = 1.0;
                        }
                    }
                }
            }
        }
        save_output_nifti(fin, "centroid_coord", centroid, false);
    }  // Centroid option is closed

    if (thinn_option == 0) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) != 0) {
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = (int)((double) *(nii_layers_data + nxy * iz + nx * ix + iy) * (double) (Nlayer_real - 1) / ((double)(Nlayer)) + 1.5);
                    }
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = Nlayer_real;
                    }
                }
            }
        }
    }
    if (thinn_option == 1) {
        Nlayer_real = Nlayer_real + 2;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) != 0) {
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = (int)((double) *(nii_layers_data + nxy * iz + nx * ix + iy) * (double) (Nlayer_real-1) / ((double) (Nlayer)) + 1.5);
                    }
                }
            }
        }
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (*(nii_layers_data + nxy * iz + nx * ix + iy) >= 2 && *(nii_layers_data + nxy * iz + nx * ix + iy)< Nlayer_real) {
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = *(nii_layers_data + nxy * iz + nx * ix + iy) - 1;
                    }
                    if (*(nii_layers_data + nxy * iz + nx * ix + iy) == Nlayer_real) {
                        *(nii_layers_data + nxy * iz + nx * ix + iy) = *(nii_layers_data + nxy * iz + nx * ix + iy) - 2;
                    }
                }
            }
        }
    }
    save_output_nifti(fin, "layers", nii_layers, true);

    cout << "  Finished." << endl;
    return 0;
}
