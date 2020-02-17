
#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_EXTREMETR: To find extreme the TR time point.\n"
    "\n"
    "    This program tries to find the maximal/minimal value of a time \n"
    "    series and writes out the timepoint of that TR vor every voxel.\n"
    "\n"
    "Usage:\n"
    "    LN_EXTREMETR -file file.nii \n"
    "\n"
    "Options:\n"
    "    -help : Show this help.\n"
    "    -file : Input time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    // nifti_image * nim_input=NULL;
    char *fin_1 = NULL, *fin_2 = NULL;
    int ac, disp_float_eg = 0;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file1\n");
                return 1;
            }
            fin_1 = argv[ac];  // Assign pointer, no string copy
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-file'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_file_1i = nifti_image_read(fin_1, 1);
    if (!nim_file_1i) {
      fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
      return 2;
    }

    // Get dimensions of input
    int sizeSlice = nim_file_1i->nz;
    int sizePhase = nim_file_1i->nx;
    int sizeRead = nim_file_1i->ny;
    int nrep =  nim_file_1i->nt;
    int nx =  nim_file_1i->nx;
    int nxy = nim_file_1i->nx * nim_file_1i->ny;
    int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

    cout << sizeSlice << " Slices | " << sizePhase << " PhaseSteps | " << sizeRead << " Read steps | " << nrep << " Timesteps " << endl;

    nifti_image * nim_file_1 = nifti_copy_nim_info(nim_file_1i);
    nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_1->nbyper = sizeof(float);
    nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
    float *nim_file_1_data = (float *) nim_file_1->data;

    // if (!fout) {
    //     fprintf(stderr, "-- no output requested \n");
    //     return 0;
    // }
    // // assign nifti_image fname/iname pair, based on output filename
    // // (request to 'check' image and 'set_byte_order' here)
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
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    nifti_image * max_file = nifti_copy_nim_info(nim_file_1);
    max_file->nt = 1;
    max_file->nvox = nim_file_1->nvox / nrep;
    max_file->datatype = NIFTI_TYPE_INT16;
    max_file->nbyper = sizeof(float);
    max_file->data = calloc(max_file->nvox, max_file->nbyper);
    short *max_file_data = (short *) max_file->data;

    nifti_image * min_file = nifti_copy_nim_info(nim_file_1);
    min_file->nt = 1;
    min_file->nvox = nim_file_1->nvox / nrep;
    min_file->datatype = NIFTI_TYPE_INT16;
    min_file->nbyper = sizeof(float);
    min_file->data = calloc(min_file->nvox, min_file->nbyper);
    short *min_file_data = (short *) min_file->data;

    double max_val = 0;
    double min_val = 1000000000;
    int TR_max = 0;
    int TR_min = 0;

    for (int islice = 0; islice < sizeSlice; ++islice) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                max_val = 0;
                min_val = 1000000000;
                for (int it = 0; it < nrep; ++it) {
                    if (*(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) > max_val) {
                        max_val = *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy);
                        TR_max = it;
                    }
                }
                for (int it = 0; it < nrep; ++it) {
                    if (*(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) < min_val) {
                        min_val = *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy);
                        TR_min = it;
                    }
                }
                *(min_file_data + nxy*islice + nx*ix + iy) =  TR_min;
                *(max_file_data + nxy*islice + nx*ix + iy) =  TR_max;
           }
        }
    }

    cout << "  Runing also until here 5... " << endl;

    string prefix_1 = "max_TR_";
    string prefix_2 = "min_TR_";
    string filename_1 = (string) (fin_1);
    string outfilename_1 = prefix_1+filename_1;
    string outfilename_2 = prefix_2+filename_1;

    cout << "  Writing as = " << outfilename_1.c_str() << " and "<< outfilename_2.c_str() << endl;

    const char *fout_1 = outfilename_1.c_str();
    if (nifti_set_filenames(max_file, fout_1 , 1, 1)) {
        return 1;
    }
    nifti_image_write(max_file);

    const char *fout_2 = outfilename_2.c_str();
    if (nifti_set_filenames(min_file, fout_2 , 1, 1)) {
        return 1;
    }
    nifti_image_write(min_file);

    // const char *fout_5 = "debug_ing.nii";
    // if (nifti_set_filenames(growfromWM0, fout_5 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(growfromWM0);

    // const char *fout_6 = "kootrGM.nii";
    // if (nifti_set_filenames(GMkoord2, fout_6 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(GMkoord2);

    // koord.autowrite("koordinaten.nii", wopts, &prot);
    cout << "  Finished." << endl;
    return 0;
}
