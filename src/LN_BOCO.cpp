
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_BOCO: This program does BOLD correction in SS-SI VASO. It does \n"
    "         the division of nulled and not nulled imaged. \n"
    "\n"
    "Usage:\n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii \n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -shift \n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 24 \n"
    "    \n"
    "    test application in the test_data folder would be:\n"
    "    ../LN_BOCO -Nulled lo_Nulled_intemp.nii -BOLD lo_BOLD_intemp.nii -trialBOCO 40 -shift\n"
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
    "    -output    : (Optional) Custom output name. \n"
    "                 including the path, if you want to write it as specific locations \n"
    "                 including the file extension: nii or nii.gz \n"
    "                 This will overwrite excisting files with the same name \n"
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
    bool use_outpath = false;
    char *fin_1 = NULL, *fin_2 = NULL, *fout = NULL;
    int ac, shift = 0;
    int trialdur = 0;
    if (argc < 2) return show_help();

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
            fout = fin_1;
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
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
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

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_2);
        return 2;
    }

    log_welcome("LN_BOCO");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int size_time = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;
    const int nr_voxels = size_time * size_z * size_y * size_x;

    // ========================================================================
    // Fix datatype issues
    nifti_image *nii_nulled = copy_nifti_as_float32(nii1);
    float *nii_nulled_data = static_cast<float*>(nii_nulled->data);
    nifti_image *nii_bold = copy_nifti_as_float32(nii2);
    float *nii_bold_data = static_cast<float*>(nii_bold->data);

    // Allocate new nifti
    nifti_image *nii_boco_vaso = copy_nifti_as_float32(nii1);
    float *nii_boco_vaso_data = static_cast<float*>(nii_boco_vaso->data);

    // ========================================================================
    // AVERAGE across Trials
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_boco_vaso_data + i) = *(nii_nulled_data + i)
                                    / (*(nii_bold_data + i));
    }

    // Clean VASO values that are unrealistic
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_boco_vaso_data + i) <= 0) {
            *(nii_boco_vaso_data + i) = 0;
        }
        if (*(nii_boco_vaso_data + i) >= 5) {
            *(nii_boco_vaso_data + i) = 5;
        }
    }

    // ========================================================================
    // Shift
    if (shift == 1) {
        nifti_image* correl_file  = nifti_copy_nim_info(nii_nulled);
        correl_file->nt = 7;
        correl_file->nvox = nii_nulled->nvox / size_time *7;
        correl_file->datatype = NIFTI_TYPE_FLOAT32;
        correl_file->nbyper = sizeof(float);
        correl_file->data = calloc(correl_file->nvox, correl_file->nbyper);
        float* correl_file_data = static_cast<float*>(correl_file->data);

        double vec_file1[size_time];
        double vec_file2[size_time];

        for (int shift = -3; shift <= 3; ++shift) {
            cout << "  Calculating shift = " << shift << endl;
            for (int j = 0; j != size_z * size_y * size_x; ++j) {
                for (int t = 3; t < size_time-3; ++t) {
                    *(nii_boco_vaso_data + nxyz * t + j) =
                        *(nii_nulled_data + nxyz * t + j)
                        / *(nii_bold_data + nxyz * (t + shift) + j);
                }
                for (int t = 0; t < size_time; ++t) {
                    vec_file1[t] = *(nii_boco_vaso_data + nxyz * t + j);
                    vec_file2[t] = *(nii_bold_data + nxyz * t + j);
                }
                *(correl_file_data + nxyz * (shift + 3) + j) =
                    ren_correl(vec_file1, vec_file2, size_time);
            }
        }

        // Get back to default
        for (int i = 0; i != nr_voxels; ++i) {
            *(nii_boco_vaso_data + i) = *(nii_nulled_data + i)
                                        / *(nii_bold_data + i);
        }

        // Clean VASO values that are unrealistic
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(nii_boco_vaso_data + i) <= 0) {
                *(nii_boco_vaso_data + i) = 0;
            }
            if (*(nii_boco_vaso_data + i) >= 2) {
                *(nii_boco_vaso_data + i) = 2;
            }
        }
        save_output_nifti(fout, "correlated", correl_file, false);
    }

    // ========================================================================
    // Trial average
    if (trialdur != 0) {
        cout << "  Doing BOLD correction after trial average..." << endl;
        cout << "    Trial duration is " << trialdur
        << ". This means there are " << (float)size_time / (float)trialdur
        <<  " trials recorded here." << endl;

        int nr_trials = size_time / trialdur;
        // Trial averave file
        nifti_image *nii_avg1 = nifti_copy_nim_info(nii1);
        nii_avg1->nt = trialdur;
        nii_avg1->nvox = nii1->nvox / size_time * trialdur;
        nii_avg1->datatype = NIFTI_TYPE_FLOAT32;
        nii_avg1->nbyper = sizeof(float);
        nii_avg1->data = calloc(nii_avg1->nvox, nii_avg1->nbyper);
        float  *nii_avg1_data  = static_cast<float*>(nii_avg1->data);

        nifti_image *nii_avg2 = nifti_copy_nim_info(nii1);
        nii_avg2->nt = trialdur;
        nii_avg2->nvox = nii1->nvox / size_time * trialdur;
        nii_avg2->datatype = NIFTI_TYPE_FLOAT32;
        nii_avg2->nbyper = sizeof(float);
        nii_avg2->data = calloc(nii_avg2->nvox, nii_avg2->nbyper);
        float  *nii_avg1_B_data  = static_cast<float*>(nii_avg2->data);

        float avg_Nulled[trialdur];
        float avg_BOLD[trialdur];

        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    for (int it = 0; it < trialdur; ++it) {
                        avg_Nulled[it] = 0;
                        avg_BOLD[it] = 0;
                    }
                    for (int it = 0; it < trialdur * nr_trials; ++it) {
                        int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                        avg_Nulled[it % trialdur] +=
                            *(nii_nulled_data + voxel_i) / nr_trials;
                        avg_BOLD[it % trialdur] +=
                            *(nii_bold_data + voxel_i) / nr_trials;
                    }

                    for (int it = 0; it < trialdur; ++it) {
                        int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                        *(nii_avg1_data + voxel_i) = avg_Nulled[it] / avg_BOLD[it];
                        *(nii_avg1_B_data + voxel_i) = avg_BOLD[it];

                    }
                }
            }
        }

        // Clean VASO values that are unrealistic
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    for (int it = 0; it < trialdur; ++it) {
                        int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;

                        if (*(nii_avg1_data + voxel_i) <= 0) {
                            *(nii_avg1_data + voxel_i) = 0;
                        }
                        if (*(nii_avg1_data + voxel_i) >= 2) {
                            *(nii_avg1_data + voxel_i) = 2;
                        }
                    }
                }
            }
        }
        save_output_nifti(fout, "VASO_trialAV_LN", nii_avg1, true);
        save_output_nifti(fout, "BOLD_trialAV_LN", nii_avg2, true);
    }
    if (!use_outpath) fout = fin_1;
    save_output_nifti(fout, "VASO_LN", nii_boco_vaso, true, use_outpath); 
    // the first argument is the string of the output file name  
    //       if there is no explicit output path given, this will be the file name of the main input data
    //       if there is an explicit output file name given, this wil be the user-defined name following the -output
    //       (inluding the path and including the file extension)
    // the second argument is the output file name prefix, that will be added to the above argument,
    //       this field is ignored, when the flag "use_outpath" (last argument) is selected.
    // the third argument is the pointer to the data set that is supposed to be written
    // the fourth argument states if, during the exectution of the program an the writing process should be logged
    //       this argument is optional with the default: TRUE
    // the fifth argument states if the output prefix (second argument) should be ignored or not
    //       this argument is optional the default: FALSE 

    cout << "  Finished." << endl;
    return 0;
}
