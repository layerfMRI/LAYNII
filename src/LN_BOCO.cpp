
#include "./common.h"
#include "./renzo_stat.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_BOCO: This program does BOLD correction in SS-SI VASO. It does \n"
    "         the division of nulled and not nulled imaged. \n"
    "\n"
    "Usage : LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii \n"
    "      : LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -shift \n"
    "      : LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 24 \n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -Nulled    : Nulled (VASO) time series that needs to be BOLD \n"
    "               : corrected.\n"
    "    -BOLD      : Reference BOLD time series without a VASO contrast.\n"
    "    -shift     : (Optional) Estimate the correlation of BOLD and VASO \n"
    "                 for temporal shifts.\n"
    "    -trialBOCO : First average trials and then do the BOLD correction. \n"
    "                 The parameter is the trial duration in TRs.\n"
    "\n"
    "Notes:\n"
    "    - Here it is assumed that BOLD and VASO refer to the double TR: \n"
    "        3dUpsample -overwrite -datum short -prefix Nulled_intemp.nii -n 2 -input Nulled.nii \n"
    "        3dUpsample -overwrite -datum short -prefix BOLD_intemp.nii -n 2 -input BOLD.nii \n"
    "    - Here I assume that they have the same spatiotemporal dimensions. \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    // nifti_image* nim_input=NULL;
    char * fin_1 = NULL, * fin_2 = NULL;
    int ac, disp_float_eg = 0, shift = 0;
    int trialdur = 0;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-Nulled")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Nulled\n");
                return 1;
            }
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-BOLD")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -BOLD\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-trialBOCO")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -trialBOCO\n");
                return 1;
            }
            trialdur = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-shift")) {
            shift = 1;
            cout << "Do a correlation analysis with temporal shifts."  << endl;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-Nulled'.\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-BOLD'.\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nim_file_1i = nifti_image_read(fin_1, 1);
    if (!nim_file_1i) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_1);
        return 2;
    }
    log_nifti_descriptives(nim_file_1i);

    nifti_image* nim_file_2i = nifti_image_read(fin_2, 1);
    if (!nim_file_2i) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_2);
        return 2;
    }
    log_nifti_descriptives(nim_file_2i);

    // Get dimensions of input
    int size_x = nim_file_1i->nx;  // phase
    int size_y = nim_file_1i->ny;  // read
    int size_z = nim_file_1i->nz;  // slice
    int size_t = nim_file_1i->nt;  // time
    int nx = nim_file_1i->nx;
    int nxy = nim_file_1i->nx * nim_file_1i->ny;
    int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

    nifti_image* nim_file_1 = nifti_copy_nim_info(nim_file_1i);
    nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_1->nbyper = sizeof(float);
    nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
    float* nim_file_1_data = (float*) nim_file_1->data;

    nifti_image* nim_file_2 = nifti_copy_nim_info(nim_file_1i);
    nim_file_2->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_2->nbyper = sizeof(float);
    nim_file_2->data = calloc(nim_file_2->nvox, nim_file_2->nbyper);
    float* nim_file_2_data = (float*) nim_file_2->data;

    // if (!fout) { fprintf(stderr, "-- no output requested \n"); return 0; }
    //     assign nifti_image fname/iname pair, based on output filename
    //     (request to 'check' image and 'set_byte_order' here)
    // if (nifti_set_filenames(nim_input, fout, 1, 1)) return 1;

    if (nim_file_1i->datatype == NIFTI_TYPE_FLOAT32) {
        float  *nim_file_1i_data = (float*) nim_file_1i->data;
        for (int it = 0; it < size_t; ++it) {
            for (int islice = 0; islice < size_z; ++islice) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz *it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_1i->datatype == NIFTI_TYPE_INT16) {
        short *nim_file_1i_data = (short *) nim_file_1i->data;
        for (int it = 0; it < size_t; ++it) {
            for (int islice = 0; islice < size_z; ++islice) {
                for (int iy=0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        *(nim_file_1_data + nxyz *it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    // if (nim_file_1i->datatype == NIFTI_TYPE_FLOAT32) {
    //     float* nim_file_1i_data = (float*) nim_file_1i->data;
    //     for (int it = 0; it < size_t; ++it) {
    //         for (int islice = 0; islice < size_z; ++islice) {
    //             for (int iy = 0; iy < size_x; ++iy) {
    //                 for (int ix = 0; ix < size_y; ++ix) {
    //                     *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_1i_data + nxyz * it + nxy * islice + nx * ix + iy));
    //                 }
    //             }
    //         }
    //     }
    // }
    if (nim_file_2i->datatype == NIFTI_TYPE_INT16) {
        short  *nim_file_2i_data = (short *) nim_file_2i->data;
        for (int it = 0; it < size_t; ++it) {
            for (int islice = 0; islice < size_z; ++islice) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        *(nim_file_2_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_2i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_file_2i->datatype == NIFTI_TYPE_FLOAT32) {
        float* nim_file_2i_data = (float*) nim_file_2i->data;
        for (int it = 0; it < size_t; ++it) {
            for (int islice = 0; islice < size_z; ++islice) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        *(nim_file_2_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_file_2i_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    float current_vaso = 0;
    nifti_image* boco_vaso = nifti_copy_nim_info(nim_file_1);
    boco_vaso->datatype = NIFTI_TYPE_FLOAT32;
    boco_vaso->nbyper = sizeof(float);
    boco_vaso->data = calloc(boco_vaso->nvox, boco_vaso->nbyper);
    float  *boco_vaso_data = (float*) boco_vaso->data;

    // AVERAGE across Trials
    for (int islice = 0; islice < size_z; ++islice) {
        for (int iy = 0; iy < size_x; ++iy) {
            for (int ix = 0; ix < size_y; ++ix) {
                for (int it = 0; it < size_t; ++it) {
                    *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = *(nim_file_1_data + nxyz * (it) + nxy * islice + nx * ix + iy) / (*(nim_file_2_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    // Clean VASO values that are unrealistic
    for (int islice = 0; islice < size_z; ++islice) {
        for (int iy = 0; iy < size_x; ++iy) {
            for (int ix=0; ix < size_y; ++ix) {
                for (int it = 0; it < size_t; ++it) {
                    if (*(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) <= 0) {
                        *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
                    }
                    if (*(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) >= 5) {
                        *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = 5;
                    }
                }
            }
        }
    }
    if (shift == 1) {
        nifti_image* correl_file  = nifti_copy_nim_info(nim_file_1);
        correl_file->nt = 7;
        correl_file->nvox = nim_file_1->nvox / size_t *7;
        correl_file->datatype = NIFTI_TYPE_FLOAT32;
        correl_file->nbyper = sizeof(float);
        correl_file->data = calloc(correl_file->nvox, correl_file->nbyper);
        float* correl_file_data = (float*) correl_file->data;

        double vec_file1[size_t];
        double vec_file2[size_t];

        for (int shift = -3; shift <= 3; ++shift) {
            cout << "  Calculating shift = " << shift << endl;
            for (int islice = 0; islice < size_z; ++islice) {
                for (int iy = 0; iy < size_x; ++iy) {
                    for (int ix = 0; ix < size_y; ++ix) {
                        for (int it = 3; it < size_t-3; ++it) {
                            *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = *(nim_file_1_data + nxyz * (it) + nxy * islice + nx * ix + iy) / (*(nim_file_2_data + nxyz * (it + shift) + nxy * islice + nx * ix + iy));
                        }
                        for (int it = 0; it < size_t; ++it) {
                            vec_file1[it] = *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy);
                            vec_file2[it] = *(nim_file_2_data + nxyz * it + nxy * islice + nx * ix + iy);
                        }
                        // TODO: WTF is `+ +`
                        *(correl_file_data + +nxyz*(shift+3) + nxy*islice + nx*ix + iy) = ren_correl(vec_file1, vec_file2,  size_t);
                    }
                }
            }
        }
        // Get back to default
        for (int islice = 0; islice < size_z; ++islice) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    for (int it = 0; it < size_t; ++it) {
                        *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = *(nim_file_1_data + nxyz * (it) + nxy * islice + nx * ix + iy) / (*(nim_file_2_data + nxyz*it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
        // Clean VASO values that are unrealistic
        for (int islice = 0; islice < size_z; ++islice) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    for (int it = 0; it < size_t; ++it) {
                        if (*(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) <= 0) {
                            *(boco_vaso_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
                        }
                        if (*(boco_vaso_data + nxyz*it + nxy*islice + nx*ix + iy) >= 2) {
                            *(boco_vaso_data + nxyz*it + nxy*islice + nx*ix + iy) = 2;
                        }
                    }
                }
            }
        }

        string prefix = "correlated_";
        string filename_1 = (string) (fin_1);
        string outfilename = prefix+filename_1;
        log_output(outfilename.c_str());
        // cout << "  Writing as = " << outfilename.c_str() << endl;

        const char* fout_1 = outfilename.c_str();
        if (nifti_set_filenames(correl_file, fout_1, 1, 1)) {
            return 1;
        }
        nifti_image_write(correl_file);
    }  // shift loop closed

    if (trialdur != 0) {
        cout << "  Also do the BOLD correction after the trial average." << endl;
        cout << "  Trial duration is " <<trialdur << ". This means there are " << (float)size_t/(float)trialdur <<  " trials recorded here." << endl;

        int numberofTrials = size_t/trialdur;
        // Trial averave file
        nifti_image* triav_file = nifti_copy_nim_info(nim_file_1);
        triav_file->nt = trialdur;
        triav_file->nvox = nim_file_1->nvox / size_t * trialdur;
        triav_file->datatype = NIFTI_TYPE_FLOAT32;
        triav_file->nbyper = sizeof(float);
        triav_file->data = calloc(triav_file->nvox, triav_file->nbyper);
        float* triav_file_data = (float*) triav_file->data;

        nifti_image* triav_B_file = nifti_copy_nim_info(nim_file_1);
        triav_B_file->nt = trialdur;
        triav_B_file->nvox = nim_file_1->nvox / size_t * trialdur;
        triav_B_file->datatype = NIFTI_TYPE_FLOAT32;
        triav_B_file->nbyper = sizeof(float);
        triav_B_file->data = calloc(triav_B_file->nvox, triav_B_file->nbyper);
        float* triav_B_file_data = (float*) triav_B_file->data;

        float AV_Nulled[trialdur];
        float AV_BOLD[trialdur];

        for (int islice = 0; islice < size_z; ++islice) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    for (int it = 0; it < trialdur; ++it) {
                        AV_Nulled[it] = 0;
                        AV_BOLD[it] = 0;
                    }
                    for (int it = 0; it < trialdur * numberofTrials; ++it) {
                        AV_Nulled[it%trialdur] = AV_Nulled[it%trialdur] + (*(nim_file_1_data + nxyz *(it) + nxy * islice + nx * ix + iy)) / numberofTrials;
                        AV_BOLD[it%trialdur] = AV_BOLD[it%trialdur] + (*(nim_file_2_data + nxyz *(it) + nxy * islice + nx * ix + iy)) / numberofTrials;
                    }
                    for (int it = 0; it < trialdur; ++it) {
                        *(triav_file_data + nxyz * it + nxy * islice + nx * ix + iy) = AV_Nulled[it] / AV_BOLD[it];
                        *(triav_B_file_data + nxyz * it + nxy * islice + nx * ix + iy) = AV_BOLD[it];
                    }
                }
            }
        }

        // clean VASO values that are unrealistic
        for (int islice = 0; islice < size_z; ++islice) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    for (int it = 0; it < trialdur; ++it) {
                        if (*(triav_file_data + nxyz * it + nxy * islice + nx * ix + iy) <= 0) {
                            *(triav_file_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
                        }
                        if (*(triav_file_data + nxyz * it + nxy * islice + nx * ix + iy) >= 2) {
                            *(triav_file_data + nxyz * it + nxy * islice + nx * ix + iy) = 2;
                        }
                    }
                }
            }
        }
        const char* fout_trial="VASO_trialAV_LN.nii";
        if (nifti_set_filenames(triav_file, fout_trial, 1, 1)) {
            return 1;
        }
        nifti_image_write(triav_file);

        const char* fout_trial_BOLD = "BOLD_trialAV_LN.nii";
        if (nifti_set_filenames(triav_B_file, fout_trial_BOLD, 1, 1)) {
            return 1;
        }
        nifti_image_write(triav_B_file);
    }  // Trial Average loop closed

    // cout << "  Running also until here 5... " << endl;
    const char* fout_5 = "VASO_LN.nii";
    log_output(fout_5);
    if (nifti_set_filenames(boco_vaso, fout_5, 1, 1)) {
        return 1;
    }
    nifti_image_write(boco_vaso);

// const char  *fout_6="kootrGM.nii";
// if (nifti_set_filenames(GMkoord2, fout_6 , 1, 1)) return 1;
// nifti_image_write(GMkoord2);

// koord.autowrite("koordinaten.nii", wopts, &prot);
    cout << "  Finished." << endl;
    return 0;
}
