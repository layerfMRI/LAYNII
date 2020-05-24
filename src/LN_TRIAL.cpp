

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_TRIAL: Average trials of block design fMRI experiments. Assumes \n"
    "          equi-distant blocks with identical rest and activity periods.\n"
    "\n"
    "Usage:\n"
    "    LN_TRIAL -input timeseries.nii -trialdur 12 \n"
    "\n"
    "for test in test folder: ../LN_TRIAL -input lo_BOLD_intemp.nii -trialdur 20 \n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -input     : Input time series.\n"
    "    -trial_dur : Duration of activity-rest trial in TRs.\n"
    "    -output    : (Optional) Custom output name. \n"
    "                 including the path, if you want to write it at specific locations \n"
    "                 including the file extension: nii or nii.gz \n"
    "                 This will overwrite excisting files with the same name \n"
    "\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ; 
    char *fin = NULL;
    int ac;
    int trial_dur;
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-trialdur")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -trial_dur\n");
                return 1;
            }
            trial_dur = atof(argv[ac]);
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

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset
    nifti_image *nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_TRIAL");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    const int size_x = nii_input->nx;
    const int size_y = nii_input->ny;
    const int size_z = nii_input->nz;
    const int size_time = nii_input->nt;
    const int nx = nii_input->nx;
    const int nxy = nii_input->nx * nii_input->ny;
    const int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    const int nr_trials = size_time / trial_dur;

    cout << "  Trial duration is " << trial_dur << ". This means there are "
         << nr_trials << " trials recorded here." << endl;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate trial average file
    nifti_image* nii_trials = nifti_copy_nim_info(nii);
    nii_trials->nt = trial_dur;
    nii_trials->nvox = nii->nvox / nr_trials;
    nii_trials->data = calloc(nii_trials->nvox, nii_trials->nbyper);
    float* nii_trials_data = static_cast<float*>(nii_trials->data);

    // ========================================================================

    for (int it = 0; it < (trial_dur * nr_trials); ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;
                    *(nii_trials_data + nxyz * (it % trial_dur) + voxel_i) +=
                        (*(nii_data + nxyz * it + voxel_i)) / nr_trials;
                }
            }
        }
    }
    
    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "TrialAverage", nii_trials, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
