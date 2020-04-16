

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_TEMPSMOOTH : Smooths data within the time domain. It removes high\n"
    "                frequency spikes like a low-pass filter.\n"
    "\n"
    "Usage:\n"
    "    LN_TEMPSMOOTH -input timeseries.nii -gaus 1.0 \n"
    "    LN_TEMPSMOOTH -input timeseries.nii -box 1 \n"
    "\n"
    "test usage in the test_data folder \n"
    "    ../LN_TEMPSMOOTH -input lo_BOLD_intemp.nii -box 3 \n"
    "\n"
    "An application of this program is described on this blog post: \n"
    "    https://layerfmri.com/anatomically-informed-spatial-smoothing/ \n"  
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii) file with time series data that will be \n"
    "              nii_smooth. Only the first time point is used. \n"
    "    -gaus   : Doing the smoothing with a Gaussian weight function. \n"
    "              A travelling window of averaging. Specify the value \n"
    "              of the Gaussian size (float values) in units of TR. \n"
    "    -box    : Doing the smoothing with a box-var. Specify the value \n"
    "              of the box sice (integer value). This is like a \n"
    "              running average sliding window.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char* fin = NULL;
    int ac, do_gaus = 0, do_box = 0, bFWHM_val = 0;
    float gFWHM_val = 0.0;
    if (argc  <  3) return show_help();

    // Process user options
    for (ac = 1; ac  <  argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-gaus")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -gaus\n");
                return 1;
            }
            gFWHM_val = atof(argv[ac]);
            do_gaus = 1;
        } else if (!strcmp(argv[ac], "-box")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -box\n");
                return 1;
            }
            bFWHM_val = atoi(argv[ac]);
            do_box = 1;
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
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
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }
    if (do_box + do_gaus !=1) {
        cout << "  Invalid smoothing option. Select gaus or box." << endl;
        return 2;
    }

    log_welcome("LN_TEMPSMOOTH");
    log_nifti_descriptives(nii_input);
    if (do_gaus) {
        cout << "Selected temporal smoothing: Gaussian" << endl;
    } else if (do_box) {
        cout << "Selected temporal smoothing: Box-car" << endl;
    }

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nr_voxels = size_x * size_y * size_z;
    // int nx = nii_input->nx;
    // int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;
    float dT = 1;

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocating necessary files
    nifti_image* nii_smooth = copy_nifti_as_float32(nii);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);

    nifti_image* nii_weight = nifti_copy_nim_info(nii);
    nii_weight->nt = 1;
    nii_weight->datatype = NIFTI_TYPE_FLOAT32;
    nii_weight->nbyper = sizeof(float);
    nii_weight->nvox = size_x * size_y * size_z;
    nii_weight->data = calloc(nii_weight->nvox, nii_weight->nbyper);
    float* nii_weight_data = static_cast<float*>(nii_weight->data);

    // ========================================================================
    // Smoothing loop
    // ========================================================================
    int vic;
    if (do_gaus) {
        vic = max(1., 2. * gFWHM_val / dT);  // Ignore if voxel is too far
    } else if (do_box) {
        vic = bFWHM_val;
    }
    cout << "    vic " << vic << endl;
    cout << "    FWHM_val " << gFWHM_val << endl;

    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_weight_data + i) = 0;
        if (*(nii_data + i) != 0) {
            for (int it = 0; it < size_time; ++it) {
                int j = nxyz * it + i;
                *(nii_smooth_data + j) = 0;

                if (do_gaus) {
                    float weight = 0;
                    int jt_start = max(0, it - vic);
                    int jt_stop = min(it + vic + 1, size_time);
                    for (int jt = jt_start; jt < jt_stop; ++jt) {
                        int k = nxyz * jt + i;
                        float dist = abs(it - jt);
                        float g = gaus(dist, gFWHM_val);
                        *(nii_smooth_data + j) += (*(nii_data + k) * g);
                        weight += g;
                    }
                    *(nii_smooth_data + j) /= weight;
                } else if (do_box) {
                    float weight = 0;
                    int jt_start = max(0, it - vic);
                    int jt_stop = min(it + vic + 1, size_time);
                    for (int jt = jt_start; jt < jt_stop; ++jt) {
                        int k = nxyz * jt + i;
                        *(nii_smooth_data + j) += *(nii_data + k);
                        weight += 1;
                    }
                    *(nii_smooth_data + j) /= weight;
                }
            }
        }
    }
    save_output_nifti(fin, "tempsmooth", nii_smooth, true);

    cout << "  Finished." << endl;
    return 0;
}
