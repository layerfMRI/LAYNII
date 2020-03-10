
// TODO(Faruk): Seems that variance flag in CLI is actually used as stdev.
// Need to ask Renzo about this.

#include "../dep/laynii_lib.h"

int N_rand = 1000;
double_t lower = -10;
double_t upper = 10;

double verteilung(double x);
typedef double (*Functions)(double);
Functions pFunc = verteilung;

double_t arb_pdf_num(int N_rand, double (*pFunc)(double),
                     double_t lower, double_t upper);
double adjusted_rand_numbers(double mean, double stdev, double value);

int show_help(void) {
    printf(
    "LN_NOISEME: Adds noise to image.\n"
    "\n"
    "Usage:\n"
    "    LN_NOISEME -input input_example.nii -variance 0.4445 \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -input      : Specify input dataset.\n"
    "    -variance   : Noise variance.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac;
    float variance_val;

    // Process user options
    if (argc < 2) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-variance")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            variance_val = atof(argv[ac]);
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image* nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_NOISEME");
    log_nifti_descriptives(nii_input);
    cout << "  Varience chosen to " << variance_val << endl;

    // ========================================================================
    // Allocating new nifti
    nifti_image* nii_new = copy_nifti_as_float32(nii_input);
    float* nii_new_data = static_cast<float*>(nii_new->data);
    // ========================================================================

    int nr_voxels = nii_input->nvox;
    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_new_data + i) +=
            adjusted_rand_numbers(0, variance_val,
                                  arb_pdf_num(N_rand, pFunc, lower, upper));
    }

    save_output_nifti(fin, "noised", nii_new, true);

    cout << "  Finished." << endl;
    return 0;
}

// Gauss lower = -5 , upper = 5
double verteilung(double z) {
    return exp(-z * z / (2.)) * 1. / sqrt(2. * M_PI);
}

double_t arb_pdf_num(int N_rand, double (*pFunc)(double), double_t lower,
                     double_t upper) {
    double_t binwidth = (upper - lower)/(double_t)N_rand;
    double_t integral = 0.0;
    double_t rand_num = rand() / (double_t)RAND_MAX;
    // rand_num = rand()/(double_t)RAND_MAX;
    // rand_num = rand()/(double_t)RAND_MAX;
    int i;

    for (i = 0; integral < rand_num; i++) {
        integral += pFunc(lower + (double_t)i * binwidth) * binwidth;

        if ((lower + (double_t)i * binwidth) > upper) {
            cout << "  Upper limit, maybe there should be adjusted limit "
                 << i << endl;
            return lower + (double_t) i *binwidth;
        }
    }
    return lower + (double_t)i * binwidth;
}

double adjusted_rand_numbers(double mean, double stdev, double value) {
    return value * stdev + mean;
}
