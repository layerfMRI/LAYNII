
#include "./common.h"
#include "./renzo_stat.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_SKEW: A program that calculates the skew of a time series.\n"
    "\n"
    "    This is helpful for artifact hunting (e.g. ghosting).\n"
    "\n"
    "Usage:\n"
    "    LN_SKEW -timeseries Nulled_intemp.nii  \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -timeseries : Any time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    // nifti_image * nim_input=NULL;
    char *fin_1 = NULL, *fin_2 = NULL;
    int ac, disp_float_eg = 0, shift = 0;
    int trialdur = 0;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }
    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-timeseries")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -timeseries\n");
                return 1;
            }
            fin_1 = argv[ac];  // Pointer assignment, no string copy
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-timeseries'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_file_1i = nifti_image_read(fin_1, 1);
    if (!nim_file_1i) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_1);
        return 2;
    }

    log_welcome("LN_SKEW");
    log_nifti_descriptives(nim_file_1i);

    // Get dimensions of input
    int sizeSlice = nim_file_1i->nz;
    int sizePhase = nim_file_1i->nx;
    int sizeRead = nim_file_1i->ny;
    int nrep = nim_file_1i->nt;
    int nx = nim_file_1i->nx;
    int nxy = nim_file_1i->nx * nim_file_1i->ny;
    int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

    nifti_image * nim_file_1 = nifti_copy_nim_info(nim_file_1i);
    nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_1->nbyper = sizeof(float);
    nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
    float *nim_file_1_data = (float *) nim_file_1->data;

    nifti_image * nim_file_2 = nifti_copy_nim_info(nim_file_1i);
    nim_file_2->datatype = NIFTI_TYPE_FLOAT32;
    nim_file_2->nbyper = sizeof(float);
    nim_file_2->data = calloc(nim_file_2->nvox, nim_file_2->nbyper);
    float  *nim_file_2_data = (float *) nim_file_2->data;

    // if (!fout) {
    //     fprintf(stderr, "-- no output requested \n");
    //     return 0;
    // }
    // // Assign nifti_image fname/iname pair, based on output filename
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

    nifti_image * skew_file = nifti_copy_nim_info(nim_file_1);
    skew_file->nt = 1;
    skew_file->nvox = nim_file_1->nvox / nrep;
    skew_file->datatype = NIFTI_TYPE_FLOAT32;
    skew_file->nbyper = sizeof(float);
    skew_file->data = calloc(skew_file->nvox, skew_file->nbyper);
    float *skew_file_data = (float *) skew_file->data;

    nifti_image * kurt_file = nifti_copy_nim_info(nim_file_1);
    kurt_file->nt = 1;
    kurt_file->nvox = nim_file_1->nvox / nrep;
    kurt_file->datatype = NIFTI_TYPE_FLOAT32;
    kurt_file->nbyper = sizeof(float);
    kurt_file->data = calloc(kurt_file->nvox, kurt_file->nbyper);
    float  *kurt_file_data = (float *) kurt_file->data;

    nifti_image * autoc_file = nifti_copy_nim_info(nim_file_1);
    autoc_file->nt = 1;
    autoc_file->nvox = nim_file_1->nvox / nrep;
    autoc_file->datatype = NIFTI_TYPE_FLOAT32;
    autoc_file->nbyper = sizeof(float);
    autoc_file->data = calloc(autoc_file->nvox, autoc_file->nbyper);
    float  *autoc_file_data = (float *) autoc_file->data;

    double vec_file1[nrep];
    double vec_file2[nrep];

    cout << "  Calculating skew, kurtosis and autocorrelation = " << shift << endl;

    for (int islice = 0; islice < sizeSlice; ++islice) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                for (int it = 0; it < nrep; ++it) {
                    vec_file1[it] = (double) * (nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy);
                }
                *(skew_file_data + nxyz * 0 + nxy * islice + nx * ix + iy) = ren_skew(vec_file1, nrep);  // gsl_stats_skew(vec_file1, 1, nrep);
                *(autoc_file_data + nxyz * 0 + nxy * islice + nx * ix + iy) = ren_autocor(vec_file1, nrep);  // gsl_stats_lag1_autocorrelation(vec_file1, 1, nrep);
                *(kurt_file_data + nxyz * 0 + nxy * islice + nx * ix + iy) = ren_kurt(vec_file1, nrep);  // gsl_stats_kurtosis(vec_file1, 1, nrep);
                // gsl_stats_skew (const double data[], size_t stride, size_t n)
                // gsl_stats_kurtosis (const double data[], size_t stride, size_t n)
                // gsl_stats_lag1_autocorrelation (const double data[], const size_t stride, const size_t n)
                // *(correl_file_data + + nxyz*(shift+3) + nxy*islice + nx*ix + iy) = gsl_stats_correlation(vec_file1, 1, vec_file2, 1, nrep);
            }
        }
    }

    string prefix_1 = "skew_";
    string filename_1 = (string) (fin_1);
    string outfilename_1 = prefix_1+filename_1;
    cout << "  Writing skew file as = " << outfilename_1.c_str() << endl;
    const char *fout_1 = outfilename_1.c_str();
    if (nifti_set_filenames(skew_file, fout_1, 1, 1)) return 1;
    nifti_image_write(skew_file);

    string prefix_2 = "kurt_";
    string filename_2 = (string) (fin_1);
    string outfilename_2 = prefix_2+filename_2;
    cout << "  Writing Kurtosis file as = " << outfilename_2.c_str() << endl;
    const char *fout_2 = outfilename_2.c_str();
    if (nifti_set_filenames(kurt_file, fout_2, 1, 1)) return 1;
    nifti_image_write(kurt_file);

    string prefix_3 = "autocorr_";
    string filename_3 = (string) (fin_1);
    string outfilename_3 = prefix_3+filename_3;
    cout << "  Writing Kurtosis file as = " << outfilename_3.c_str() << endl;
    const char *fout_3 = outfilename_3.c_str();
    if (nifti_set_filenames(autoc_file, fout_3, 1, 1)) return 1;
    nifti_image_write(autoc_file);

    /////////////////////////////////////////////
    // Calculating correlation with everything //
    /////////////////////////////////////////////
    nifti_image * conc_file = nifti_copy_nim_info(nim_file_1);
    conc_file->nt = 1;
    conc_file->nvox = nim_file_1->nvox / nrep;
    conc_file->datatype = NIFTI_TYPE_FLOAT32;
    conc_file->nbyper = sizeof(float);
    conc_file->data = calloc(conc_file->nvox, conc_file->nbyper);
    float *conc_file_data = (float *) conc_file->data;

    for (int it = 0; it < nrep; ++it) {
        vec_file1[it] = 0;
        vec_file2[it] = 0;
    }

    // cout << " number of voxerls"  <<nxyz << endl;
    // Mean time course of everything
    for (int islice = 0; islice < sizeSlice; ++islice) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                for (int it = 0; it < nrep; ++it) {
                    vec_file1[it] = vec_file1[it] + (double) *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) / nxyz;
                    // if ( *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy) / nxyz > 10000) cout << "  I am weird " << endl;
                }
            }
        }
    }
    // Voxel-wise corelation to mean of everything
    for (int islice = 0; islice < sizeSlice; ++islice) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix <sizeRead; ++ix) {
                for (int it = 0; it < nrep; ++it)   {
                    vec_file2[it] = (double) *(nim_file_1_data + nxyz * it + nxy * islice + nx * ix + iy);
                }
                *(conc_file_data + nxyz * 0 + nxy * islice + nx * ix + iy) = ren_correl(vec_file1, vec_file2, nrep);
            }
        }
    }

    string prefix_4 = "overall_correl_";
    string filename_4 = (string) (fin_1);
    string outfilename_4 = prefix_4+filename_4;
    cout << "  Writing Overall correlation file as = " << outfilename_4.c_str() << endl;
    const char *fout_4 = outfilename_4.c_str();
    if (nifti_set_filenames(conc_file, fout_4, 1, 1)) {
        return 1;
    }
    nifti_image_write(conc_file);

    // const char *fout_6 = "kootrGM.nii";
    // if (nifti_set_filenames(GMkoord2, fout_6, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(GMkoord2);

    // koord.autowrite("koordinaten.nii", wopts, &prot);
    cout << "  Finished." << endl;
    return 0;
}
