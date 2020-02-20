

#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_SMOOTH_RIM : Corrects errors in the manually drawn rim file by \n"
    "                filling holes and removing voxels that stick out of \n"
    "                the smooth surfaces. \n"
    "\n"
    "Usage:\n"
    "    LN_SMOOTH_RIM -rim rim.nii \n"
    "\n"
    "Options:\n"
    "   -help : Show this help.\n"
    "   -rim  : File of voxels with 1 (CSF border), 2 (WM border), 3 (GM ribbon).\n"
    "   -mask : (Optional) Mask activity outside of layers. \n"
    "\n"
    "Note: This program now supports INT16, INT32 and FLOAT23 \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *rim_filename = NULL, *fout = NULL;
    int ac, twodim = 0, do_masking = 0;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
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
            rim_filename = argv[ac]; // no string copy, just pointer assignment
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!rim_filename) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_rim_r = nifti_image_read(rim_filename, 1);
    if (!nim_rim_r) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", rim_filename);
        return 2;
    }

    log_welcome("LN_SMOOTH_RIM");
    log_nifti_descriptives(nim_rim_r);

    // Get dimensions of input
    int sizeSlice = nim_rim_r->nz;
    int sizePhase = nim_rim_r->nx;
    int sizeRead = nim_rim_r->ny;
    int nrep = nim_rim_r->nt;
    int nx = nim_rim_r->nx;
    int nxy = nim_rim_r->nx * nim_rim_r->ny;
    int nxyz = nim_rim_r->nx * nim_rim_r->ny * nim_rim_r->nz;
    float dX = nim_rim_r->pixdim[1];
    float dY = nim_rim_r->pixdim[2];
    float dZ = nim_rim_r->pixdim[3];

    if (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_rim = nifti_copy_nim_info(nim_rim_r);
    nim_rim->datatype = NIFTI_TYPE_INT16;
    nim_rim->nbyper = sizeof(short);
    nim_rim->data = calloc(nim_rim->nvox, nim_rim->nbyper);
    short  *nim_rim_data = (short *) nim_rim->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_rim_r->datatype == NIFTI_TYPE_FLOAT32) {
        float  *nim_rim_r_data = (float *) nim_rim_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_rim_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_rim_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_rim_r->datatype == NIFTI_TYPE_INT16) {
        short  *nim_rim_r_data = (short *) nim_rim_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_rim_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_rim_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_rim_r->datatype == NIFTI_TYPE_INT32) {
        int *nim_rim_r_data = (int *) nim_rim_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_rim_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_rim_r_data + nxyz * it + nxy * islice + nx * ix + iy));
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
                if (*(nim_rim_data + nxy * iz + nx * ix + iy) > layernumber) layernumber = *(nim_rim_data + nxy * iz + nx * ix + iy);
            }
        }
    }

    cout << "  There are " << layernumber << " to smooth." << endl;
    if (layernumber > 3) {
        cout << "  ################################### " << endl;
        cout << "  ##       You are in Trouble      ## " << endl;
        cout << "  ## There are more than 3 numbers ## " << endl;
        cout << "  ################################### " << endl;
    }

    //////////////////////////////////////////////
    // Removing voxels that are too deep in CSF //
    //////////////////////////////////////////////
    float vox_dist(float x1, float y1, float z1, float x2, float y2, float z2);
    int vinc = 1;
    int isedge = 0;
    float dist_i = 0;

    for (int ilayer = 1; ilayer <=layernumber; ++ilayer) {
        cout << "  ilayer " << endl;

        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    //  *(smoothed_data + nxy * iz + nx * ix + iy) = 0 ;
                    // Remove warts
                    if (*(nim_rim_data + nxy * iz + nx * ix + iy) == ilayer) {
                        isedge = 0;
                        for (int iz_i = max(0, iz - vinc); iz_i < min(iz+vinc+1, sizeRead); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc); iy_i < min(iy+vinc+1, sizePhase); ++iy_i) {
                                for (int ix_i = max(0, ix - vinc); ix_i < min(ix+vinc+1, sizeRead); ++ix_i) {
                                    dist_i = vox_dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i);
                                    if (*(nim_rim_data + nxy*iz_i + nx*ix_i + iy_i) == 0 &&  dist_i == 1) {
                                        isedge = isedge + 1;
                                    }
                                }
                            }
                        }
                        if (isedge > 4) {
                            *(nim_rim_data + nxy*iz + nx*ix + iy) = 0;
                        }
                    }
                }
            }
        }
    }  // Layer loop closed

    cout << "  Writing output..." << endl;
    string prefix_1 = "surfacesmoothed_";
    string filename_1 = (string) (rim_filename);
    string outfilename_1 = prefix_1+filename_1;
    cout << "  Writing as = " << outfilename_1.c_str() << endl;

    const char *fout_1 = outfilename_1.c_str();
    if (nifti_set_filenames(nim_rim, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_rim);

    // const char *fout_1 = "layer.nii" ;
    // if (nifti_set_filenames(layer, fout_1 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(layer);

    cout << "  Finished." << endl;
    return 0;
}

float vox_dist(float x1, float y1, float z1, float x2, float y2, float z2) {
    float return_value = 0;
    if (sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) <= 1.) {
        return_value = 1;
    } else {
        return_value = 0;
    }
    return return_value;
}
