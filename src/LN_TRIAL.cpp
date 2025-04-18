#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_TRIAL: Average trials of block design fMRI experiments. Assumes \n"
    "          equi-distant blocks with identical rest and activity periods.\n"
    "\n"
    "Usage:\n"
    "    LN_TRIAL -input timeseries.nii -trialdur 12 \n"
    "    ../LN_TRIAL -input lo_BOLD_intemp.nii -trialdur 20 \n" 
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -input    : Input time series.\n"
    "    -trialdur : Duration of activity-rest trial in volumes (TRs).\n"
    "    -output   : (Optional) Output filename, including .nii or .nii.gz, and\n"
    "                path if needed. Overwrites existing files.\n"    
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fout = NULL ;
    char *fin = NULL;
    int ac;
    uint64_t trial_dur;
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
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-trialdur")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -trialdur\n");
                return 1;
            }
            trial_dur = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
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
    const uint64_t size_x = nii_input->nx;
    const uint64_t size_y = nii_input->ny;
    const uint64_t size_z = nii_input->nz;
    const uint64_t size_time = nii_input->nt;

    const uint64_t nr_voxels = size_x * size_y * size_z;
    const uint64_t nr_trials = size_time / trial_dur;

    cout << "  Trial duration is " << trial_dur << ". There are " << nr_trials << " trials." << endl;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate trial average file
    nifti_image* nii_trials = nifti_copy_nim_info(nii);
    nii_trials->nt = trial_dur;
    nii_trials->nvox = nr_voxels / nr_trials;
    nii_trials->data = calloc(nii_trials->nvox, nii_trials->nbyper);
    float* nii_trials_data = static_cast<float*>(nii_trials->data);

    // ========================================================================
    for (uint64_t t = 0; t < (trial_dur * nr_trials); ++t) {
        for (uint64_t z = 0; z < size_z; ++z) {
            for (uint64_t y = 0; y < size_y; ++y) {
                for (uint64_t x = 0; x < size_x; ++x) {
                    uint64_t voxel_i = size_x * size_y * z + size_x * y + x;
                    *(nii_trials_data + nr_voxels*(t%trial_dur) + voxel_i) += 
                        (*(nii_data + nr_voxels*t + voxel_i)) / static_cast<float>(nr_trials);
                }
            }
        }
    }

    save_output_nifti(fout, "TrialAverage", nii_trials, true);

    cout << "  Finished." << endl;
    return 0;
}
