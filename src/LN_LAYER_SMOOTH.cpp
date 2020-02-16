
#include "./common.h"

int show_help(void) {
    printf(
    "LN_LAYER_SMOOTH : Layering algorithm based on iterative smothing.\n"
    "\n"
    "    This program smooths data within layer or columns. In order to \n"
    "    avoid smoothing across masks, a crawler smooths only across \n"
    "    connected voxels.\n"
    "\n"
    "Usage:\n"
    "    LN_LAYER_SMOOTH -layer_file layers.nii -input activity_map.nii -FWHM 1\n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -layer_file : Nifti (.nii) file that contains layer or column masks.\n"
    "    -input      : Nifti (.nii) file that should be smoothed. It \n"
    "                  should have same dimensions as layer file.\n"
    "    -FWHM       : The amount of smoothing in mm.\n"
    "    -twodim     : (Optional) Smooth only in 2D. \n"
    "    -mask       : (Optional) Mask activity outside of layers. \n"
    "    -sulctouch  : (Optional) Allows smoothing across sucli. This is \n"
    "                  necessary, when you do heavy smoothing well bevond \n"
    "                  the spatial scale of the cortical thickness, or heavy\n"
    "                  curvature. It will make things slower. Note that this \n"
    "                  is best done with not too many layers. Otherwise a \n"
    "                  single layer has holes and is not connected.\n"
    "\n"
    "Note: This program now supports INT16, INT32 and FLOAT32.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * fmaski = NULL, *fout = NULL, * finfi = NULL;
    int ac, twodim = 0, do_masking = 0, sulctouch = 0;
    float FWHM_val = 0;
    if (argc < 3) return show_help();  // Typing '-help' is sooo much work

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layer_file\n");
                return 1;
            }
            fmaski = argv[ac];  // No string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);  // No string copy, pointer assignment
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            finfi = argv[ac];  // No string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << "Smooth only in 2D."  << endl;
        } else if (!strcmp(argv[ac], "-sulctouch")) {
            sulctouch = 1;
            cout << "Smooth across sluci, might take longer."  << endl;
        } else if (!strcmp(argv[ac], "-mask")) {
            do_masking = 1;
            cout << "Set voxels to 0 outside layers (masking option)"  << endl;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finfi) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    // read input dataset, including data
    nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
    if (!nim_inputfi) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", finfi);
        return 2;
    }
    if (!fmaski) {
        fprintf(stderr, "** missing option '-layer_file'\n");
        return 1;
    }
    // read input dataset, including data
    nifti_image * nim_maski = nifti_image_read(fmaski, 1);
    if (!nim_maski) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", fmaski);
        return 2;
    }

    // Get dimensions of input
    int sizeSlice = nim_maski->nz;
    int sizePhase = nim_maski->nx;
    int sizeRead = nim_maski->ny;
    int nrep = nim_maski->nt;
    int nx = nim_maski->nx;
    int nxy = nim_maski->nx * nim_maski->ny;
    int nxyz = nim_maski->nx * nim_maski->ny * nim_maski->nz;
    float dX = nim_maski->pixdim[1];
    float dY = nim_maski->pixdim[2];
    float dZ = nim_maski->pixdim[3];

    if (twodim == 1) dZ = 1000 * dZ;
    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_inputf = nifti_copy_nim_info(nim_inputfi);
    nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
    nim_inputf->nbyper = sizeof(float);
    nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
    float* nim_inputf_data = (float*) nim_inputf->data;

    nifti_image * nim_mask = nifti_copy_nim_info(nim_maski);
    nim_mask->datatype = NIFTI_TYPE_INT32;
    nim_mask->nbyper = sizeof(int);
    nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
    int* nim_mask_data = (int*) nim_mask->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////

    if (nim_inputfi->datatype == NIFTI_TYPE_FLOAT32) {
        float* nim_inputfi_data = (float*) nim_inputfi->data;
        for (int it=0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_inputfi->datatype == NIFTI_TYPE_INT16) {
        short *nim_inputfi_data = (short *) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix=0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_inputfi->datatype == NIFTI_TYPE_INT32) {
        int* nim_inputfi_data = (int*) nim_inputfi->data;
        for (int it=0; it < nrep; ++it) {
            for (int islice=0; islice < sizeSlice; ++islice) {
                for (int iy=0; iy < sizePhase; ++iy) {
                    for (int ix=0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_maski->datatype == NIFTI_TYPE_FLOAT32) {
        float* nim_maski_data = (float*) nim_maski->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_mask_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_maski_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_maski->datatype == NIFTI_TYPE_INT16) {
        short *nim_maski_data = (short *) nim_maski->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_mask_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_maski_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_maski->datatype == NIFTI_TYPE_INT32) {
        int* nim_maski_data = (int*) nim_maski->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_mask_data + nxyz * it + nxy * islice + nx * ix + iy) = (int) (*(nim_maski_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    cout << sizeSlice << " Slices | " << sizePhase << " Phase steps | " << sizeRead << " Read_steps | " << nrep << " Time_steps " << endl;
    cout << "  Voxel size = " << dX << " x " << dY << " x "  << dZ << endl;
    cout << "  Datatype 1 + " << nim_inputf->datatype << endl;
    cout << "  Datatype 2 + " << nim_mask->datatype << endl;

    /////////////////////////////////////
    // MAKE allocating necessary files //
    /////////////////////////////////////

    nifti_image * smoothed = nifti_copy_nim_info(nim_inputf);
    nifti_image * gausweight = nifti_copy_nim_info(nim_inputf);

    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;

    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);

    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);

    float* smoothed_data = (float*) smoothed->data;
    float* gausweight_data = (float*) gausweight->data;

    float dist(float x1, float y1, float z1, float x2, float y2, float z2, float dX, float dY, float dZ);
    float gaus(float distance, float sigma);

    cout << "  Debug 2" << endl;

    // float kernel_size = 10;  // Corresponds to one voxel sice
    int vinc = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
    float dist_i = 0.;
    cout << "  vinc " <<  vinc <<  endl;
    cout << "  FWHM_val " <<  FWHM_val <<  endl;

    //////////////////////////////
    // Finding number of layers //
    //////////////////////////////
    int layernumber = 0;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_mask_data + nxy * iz + nx * ix + iy) > layernumber) layernumber = *(nim_mask_data + nxy * iz + nx * ix + iy);
            }
        }
    }

    cout << "  There are  " << layernumber << " layers/masks to smooth within." << endl;

    ///////////////////////
    // Making it integer //
    ///////////////////////

    ////////////////////
    // SMOOTHING LOOP //
    ////////////////////
    // cout << " DEBUG " <<   dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl;

    if (sulctouch == 0) {
        cout << "  Smoothing in layer not considering sulci." << endl;

        for (int layernumber_i = 1; layernumber_i <= layernumber; ++layernumber_i) {
            cout << "\r    " << layernumber_i << " of " << layernumber << flush;

            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead - 0; ++ix) {
                        *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                        // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;

                        if (*(nim_mask_data + nxy*iz + nx * ix + iy) == layernumber_i) {

                            for (int iz_i = max(0, iz-vinc); iz_i < min(iz + vinc + 1, sizeSlice-1); ++iz_i) {
                                for (int iy_i = max(0, iy-vinc); iy_i < min(iy + vinc + 1, sizePhase - 1); ++iy_i) {
                                    for (int ix_i = max(0, ix-vinc); ix_i < min(ix + vinc + 1, sizeRead - 1); ++ix_i) {
                                        if (*(nim_mask_data + nxy * iz_i + nx * ix_i + iy_i) == layernumber_i) {
                                            dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                            // cout << "debug  4 " << gaus(dist_i ,FWHM_val) <<   endl;
                                            // cout << "debug  5 " << dist_i  <<   endl;
                                            // if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) cout << "debug  4b " << endl;
                                            // dummy = *(layer_data  + nxy*iz_i + nx*ix_i  + iy_i );
                                            *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                            *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                        }
                                    }
                                }
                            }
                            if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) / *(gausweight_data + nxy * iz + nx * ix + iy);
                        }
                        if (*(nim_mask_data + nxy * iz + nx * ix + iy) <= 0) *(smoothed_data + nxy * iz + nx * ix + iy) = *(nim_inputf_data + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }  // for layer loop closed
        cout << endl;
    }  // if closed
    //////////////////////////////////////////////////////////////////
    // if requested I do the smoothing only within connected layers //
    //////////////////////////////////////////////////////////////////

    if (sulctouch == 1) {
        // Allocating local connected vincinity file
        nifti_image * hairy_brain = nifti_copy_nim_info(nim_mask);
        hairy_brain->datatype = NIFTI_TYPE_INT32;
        hairy_brain->nbyper = sizeof(int);
        hairy_brain->data = calloc(hairy_brain->nvox, hairy_brain->nbyper);
        int* hairy_brain_data = (int*) hairy_brain->data;
        hairy_brain->scl_slope = 1.;
        int vinc_steps = 1;

        // Making sure I am cooking in a clean kitchen
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead - 0; ++ix) {
                    *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                }
            }
        }

        vinc = max(1., 2. * FWHM_val / dX);  // Ignore if voxel is too far
        int layernumber_i = 0;  // Running index
        cout << " vinc " << vinc << endl;
        cout << " FWHM_val " << FWHM_val << endl;
        cout << " starting within sulucal smoothing now  " <<  endl;

        // For estimation of time
        int nvoxels_to_go_across = 0;
        int running_index = 0;
        int pref_ratio = 0;

        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (*(nim_mask_data + nxy * iz + nx * ix + iy) > 1) nvoxels_to_go_across++;
                }
            }
        }
        // cout << "  here 1" << endl;

        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead - 0; ++ix) {
                    if (*(nim_mask_data + nxy * iz + nx * ix + iy) > 0) {
                        running_index++;
                        if ((running_index * 100) / nvoxels_to_go_across != pref_ratio) {
                            cout << "\r " << (running_index * 100) / nvoxels_to_go_across <<  "% " << flush;
                            pref_ratio = (running_index * 100) / nvoxels_to_go_across;
                        }
                        layernumber_i = *(nim_mask_data + nxy * iz + nx * ix + iy);
                        *(gausweight_data + nxy * iz + nx * ix + iy) = 0;

                        /////////////////////////////////////////////////
                        // Find area that is not from the other sulcus //
                        /////////////////////////////////////////////////

                        // PREPARATION OF DUMMY VINSINITY FILE
                        // cout << "  here 2" << endl;

                        for (int iz_i = max(0, iz - vinc - vinc_steps); iz_i <= min(iz + vinc + vinc_steps, sizeSlice - 1); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc - vinc_steps); iy_i <= min(iy + vinc + vinc_steps, sizePhase - 1); ++iy_i) {
                                for (int ix_i = max(0, ix - vinc - vinc_steps); ix_i <= min(ix + vinc + vinc_steps, sizeRead - 1); ++ix_i) {
                                    // cout << iz_i << " " << iy_i << "  " << ix_i << "  " <<  sizeSlice-1 << " " << sizePhase-1 << "  " << sizePhase-1 << "  " << endl;
                                    *(hairy_brain_data  + nxy*iz_i + nx*ix_i  + iy_i) = 0;
                                }
                            }
                        }
                        // cout << "  here " << endl;
                        *(hairy_brain_data  + nxy*iz + nx*ix  + iy) = 1;

                        // Growing into neigbouring voxels.
                        for (int K_= 0; K_< vinc; K_++) {
                            for (int iz_ii = max(0, iz - vinc); iz_ii <= min(iz + vinc, sizeSlice - 1); ++iz_ii) {
                                for (int iy_ii = max(0,iy - vinc); iy_ii <= min(iy + vinc, sizePhase - 1); ++iy_ii) {
                                    for (int ix_ii = max(0, ix-vinc); ix_ii <= min(ix + vinc, sizeRead - 1); ++ix_ii) {
                                        if (*(hairy_brain_data + nxy * iz_ii + nx * ix_ii + iy_ii) == 1) {
                                            for (int iz_i = max(0, iz_ii - vinc_steps); iz_i <= min(iz_ii + vinc_steps, sizeSlice - 1); ++iz_i) {
                                                for (int iy_i = max(0, iy_ii - vinc_steps); iy_i<= min(iy_ii + vinc_steps, sizePhase - 1); ++iy_i) {
                                                    for (int ix_i = max(0, ix_ii - vinc_steps); ix_i <= min(ix_ii + vinc_steps, sizeRead - 1); ++ix_i) {
                                                        if (dist((float)ix_ii, (float)iy_ii, (float)iz_ii, (float)ix_i, (float)iy_i, (float)iz_i, 1, 1, 1) <= 1.74 && *(nim_mask_data + nxy * iz_i + nx * ix_i + iy_i) == layernumber_i) {
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

                        // Apply smoothing within each layer and within local patch
                        for (int iz_i = max(0, iz-vinc); iz_i <= min(iz + vinc, sizeSlice - 1); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc); iy_i <= min(iy + vinc, sizePhase - 1); ++iy_i) {
                                for (int ix_i = max(0, ix-vinc); ix_i <= min(ix + vinc, sizeRead - 1); ++ix_i) {
                                    if (*(hairy_brain_data + nxy * iz_i + nx * ix_i + iy_i) == 1) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, FWHM_val);
                                        *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                    }
                                }
                            }
                        }
                        if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) / *(gausweight_data + nxy * iz + nx * ix + iy);
                    } // if scope  if (*(nim_mask_data + nxy * iz + nx * ix + iy) > 0){ closed

                }
            }
        }
        // for(int iz=0; iz<sizeSlice; ++iz){
        //     for(int iy=0; iy<sizePhase; ++iy){
        //         for(int ix=0; ix<sizeRead; ++ix){
        //             if(*(nim_inputf_data + nxy * iz + nx * ix + iy) > 0) {
        //                 *(nim_inputf_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy);
        //             }
        //         }
        //     }
        // }
        const char *fout_3 = "hairy_brain.nii";
        if (nifti_set_filenames(hairy_brain, fout_3, 1, 1)) return 1;
        nifti_image_write(hairy_brain);
    }
    cout << "  Smoothing done. " <<  endl;

    ////////////////////////////////
    // Masking if it is it wanted //
    ////////////////////////////////

    if (do_masking == 1) {
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        if (*(nim_mask_data + nxyz * it + nxy * islice + nx * ix + iy) == 0) {
                            *(smoothed_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
                        }
                    }
                }
            }
        }
    }

    // cout << " runing also until here  5.... " << endl;
    // cout << " slope " << smoothed->scl_slope << " " << nim_inputfi->scl_slope << endl;

    smoothed->scl_slope = nim_inputfi->scl_slope;

    if (nim_inputfi->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    // Output file name
    // const char *fout_4 = "leaky_layers.nii" ;
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) return 1;
    // nifti_image_write(leak_layer);

    // const char *fout_5 = "input_file.nii" ;
    // if (nifti_set_filenames(nim_inputf, fout_5, 1, 1)) return 1;
    // nifti_image_write(nim_inputf);

    // const char *fout_2 = "mask.nii" ;
    // if (nifti_set_filenames(nim_mask, fout_2, 1, 1)) return 1;
    // nifti_image_write(nim_mask);

    string prefix = "smoothed_";
    string filename = (string) (finfi);
    string outfilename = prefix + filename;

    cout << "  Writing as = " << outfilename.c_str() << endl;  // finfi is: char *

    const char *fout_1 = outfilename.c_str();
    if (nifti_set_filenames(smoothed, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(smoothed);

    // const char  *fout_1="layer.nii" ;
    // if(nifti_set_filenames(layer, fout_1 , 1, 1)) return 1;
    // nifti_image_write(layer);
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt((x1 - x2) * (x1 - x2) * dX * dX + (y1 - y2) * (y1 - y2)
                * dY * dY + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592))
           * exp(-0.5 * distance * distance / (sigma * sigma));
}
