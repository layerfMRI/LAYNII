
// TODO(@Faruk): Compiles but gives error:
// `** ERROR: nifti_image_write_hdr_img: no image data`
// Also seems to give the same error in the old version. Check with Renzo.

#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_ZOOM: Reduces the nii file to only include the rectangular space \n"
    "         that is different from zero.\n"
    "\n"
    "Usage: \n"
    "    LN_ZOOM -mask mask.nii -input file_to_be_zoomed.nii \n"
    "\n"
    "Options:\n"
    "\n"
    "    -help  : Show this help.\n"
    "    -mask  : Nifti (.nii) file that determines the region of interest\n"
    "             (e.g. the layer mask with one time point).\n"
    "    -input : Nifti (.nii) file that should be zoomed (e.g. with \n"
    "             multiple time points).\n"
    "\n"
    "Note: Written for Insub. \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    // nifti_image * nim_input=NULL;
    char *fin_1 = NULL, *fin_2 = NULL;
    int ac, disp_float_eg = 0, shift = 0;
    int trialdur = 0;
    if (argc < 2) return show_help();   // Typing '-help' is sooo much work

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Nulled\n");
            return 1;
            }
            fin_1 = argv[ac];  // no string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -BOLD\n");
                return 1;
            }
            fin_2 = argv[ac];  // no string copy, just pointer assignment
        } else {
            fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!fin_1) {
       fprintf(stderr, "** missing option '-input'\n");
       return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_file_1i = nifti_image_read(fin_1, 1);
    if (!nim_file_1i) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
        return 2;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-mask'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_file_2i = nifti_image_read(fin_2, 1);
    if (!nim_file_2i) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_2);
        return 2;
    }
    // Get dimensions of input
    int sizeSlice = nim_file_1i->nz;
    int sizePhase = nim_file_1i->nx;
    int sizeRead = nim_file_1i->ny;
    int nrep = nim_file_1i->nt;
    int nx = nim_file_1i->nx;
    int nxy = nim_file_1i->nx * nim_file_1i->ny;
    int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

    cout << sizeSlice << " Slices " << sizePhase << " | PhaseSteps " << sizeRead << " | Read steps " << nrep << " timesteps " << endl;

    nifti_image * nim_file_1 = nifti_copy_nim_info(nim_file_1i);
    nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_1->nbyper = sizeof(float);
    nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
    float  *nim_file_1_data = (float *) nim_file_1->data;

    nifti_image * nim_file_2 = nifti_copy_nim_info(nim_file_1i);
    nim_file_2->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_2->nbyper = sizeof(float);
    nim_file_2->data = calloc(nim_file_2->nvox, nim_file_2->nbyper);
    float  *nim_file_2_data = (float *) nim_file_2->data;

    // if (!fout) {
    //     fprintf(stderr, "-- no output requested \n");
    //     return 0;
    // }
    // assign nifti_image fname/iname pair, based on output filename
    // (request to 'check' image and 'set_byte_order' here)
    // if (nifti_set_filenames(nim_input, fout, 1, 1)) {
    //     return 1;
    // }

    if (nim_file_1i->datatype == NIFTI_TYPE_FLOAT32) {
        float *nim_file_1i_data = (float *) nim_file_1i->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_1i->datatype == NIFTI_TYPE_INT16) {
        short *nim_file_1i_data = (short *) nim_file_1i->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_1i->datatype == NIFTI_TYPE_FLOAT32) {
        float *nim_file_1i_data = (float *) nim_file_1i->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy <sizePhase; ++iy) {
                    for (int ix = 0;  ix < sizeRead; ++ix) {
                        *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_2i->datatype == NIFTI_TYPE_INT16) {
        short  *nim_file_2i_data = (short *) nim_file_2i->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_file_2_data + nxyz * 0 + nxy * islice + nx * ix + iy) = (float) (*(nim_file_2i_data + nxyz * 0 + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_2i->datatype == NIFTI_TYPE_FLOAT32) {
        float *nim_file_2i_data = (float *) nim_file_2i->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_file_2_data + nxyz * 0 + nxy * islice + nx * ix + iy) = (float) (*(nim_file_2i_data + nxyz * 0 + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    int min_z = 10000;
    int max_z = 0;
    int min_x = 10000;
    int max_x = 0;
    int min_y = 10000;
    int max_y = 0;
    float mask_val = 0;

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                mask_val = *(nim_file_2_data + nxyz *0 + nxy*iz + nx*ix + iy);
                if (mask_val > 0 && ix > max_x) max_x = ix;
                if (mask_val > 0 && iy > max_y) max_y = iy;
                if (mask_val > 0 && iz > max_z) max_z = iz;
                if (mask_val > 0 && ix < min_x) min_x = ix;
                if (mask_val > 0 && iy < min_y) min_y = iy;
                if (mask_val > 0 && iz < min_z) min_z = iz;
            }
        }
    }
    cout << "  x range is " << min_x << "-" << max_x << endl;
    cout << "  y range is " << min_y << "-" << max_y << endl;
    cout << "  z range is " << min_z << "-" << max_z << endl;
    int zoomed_z_size = max_z - min_z;
    int zoomed_x_size = max_x - min_x;
    int zoomed_y_size = max_y - min_y;

    nifti_image * zoomed_file = nifti_copy_nim_info(nim_file_1);
    // zoomed_file->nt = 7;
    zoomed_file->nz = zoomed_z_size;
    zoomed_file->nx = zoomed_x_size;
    zoomed_file->ny = zoomed_y_size;

    zoomed_file->nvox = (nim_file_1->nvox * zoomed_x_size  * zoomed_y_size * zoomed_z_size) / (nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz);
    zoomed_file->datatype = NIFTI_TYPE_FLOAT32;
    zoomed_file->nbyper = sizeof(float);
    zoomed_file->data = calloc(zoomed_file->nvox, zoomed_file->nbyper);
    float *zoomed_file_data = (float *) zoomed_file->data;

    cout << "  Reduction " << (nim_file_1->nvox * zoomed_x_size  * zoomed_y_size * zoomed_z_size) / (nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz) << endl;

    int nx_z = zoomed_file->nx;
    int nxy_z = zoomed_file->nx * zoomed_file->ny;
    int nxyz_z = zoomed_file->nx * zoomed_file->ny * zoomed_file->nz;

    for (int it = 0; it < nrep; ++it) {
        for (int iz = min_z; iz < max_z; ++iz) {
            for (int iy = min_y; iy < max_y; ++iy) {
                for (int ix = min_x; ix < max_x; ++ix) {
                    *(zoomed_file_data + nxyz_z * it + nxy_z * (iz-min_z) + nx_z * (ix - min_x) + (iy - min_y)) = *(nim_file_1_data + nxyz * it + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    string prefix = "zoomed_";
    string filename_1 = (string) (fin_1);
    string outfilename = prefix+filename_1;
    cout << "  Writing as = " << outfilename.c_str() << endl;

    const char *fout_1 = outfilename.c_str();
    if (nifti_set_filenames(zoomed_file, fout_1 , 1, 1)) {
        return 1;
    }
    nifti_image_write(zoomed_file);
    cout << "  Running until here... " << endl;

    // const char  *fout_6="kootrGM.nii";
    // if (nifti_set_filenames(GMkoord2, fout_6 , 1, 1)) return 1;
    // nifti_image_write(GMkoord2);

    // koord.autowrite("koordinaten.nii", wopts, &prot);
    return 0;
}
