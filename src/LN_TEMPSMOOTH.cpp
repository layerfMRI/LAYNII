

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
            cout << "Selected temporal smoothing: Gaussian" << endl;
        } else if (!strcmp(argv[ac], "-box")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -box\n");
                return 1;
            }
            bFWHM_val = atoi(argv[ac]);
            do_box = 1;
            cout << "Selected temporal smoothing: Box-car" << endl;
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

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_t = nii_input->nt;
    int nr_voxels = size_t * size_z * size_y * size_x;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;
    float dX = nii_input->pixdim[1];
    // float dY = nii_input->pixdim[2];
    // float dZ = nii_input->pixdim[3];

    // Note: If you are running the smoothing in 2D, it will still go through
    // the entire pipeline. The only difference is that the weights in a
    // certain direction are suppressed. Doing it in 2D, will be faster.

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocating necessary files
    nifti_image* nii_smooth = copy_nifti_as_float32(nii);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);

    nifti_image* nii_gaussw = nifti_copy_nim_info(nii);
    nii_gaussw->nt = 1;
    nii_gaussw->datatype = NIFTI_TYPE_FLOAT32;
    nii_gaussw->nbyper = sizeof(float);
    nii_gaussw->nvox = nii_gaussw->nvox / size_t;
    nii_gaussw->data = calloc(nii_gaussw->nvox, nii_gaussw->nbyper);
    float* nii_gaussw_data = static_cast<float*>(nii_gaussw->data);

    // ========================================================================
    // Smoothing loop (for Gaussian smoothing
    // ========================================================================
    if (do_gaus) {
        // float_kernel_size = 10;  // Corresponds to one voxel size.
        int vinc = max(1., 2. * gFWHM_val / dX);  // Ignore if voxel is too far
        cout << "  vinc " << vinc << endl;
        cout << "  FWHM_val " << gFWHM_val << endl;

        // To estimate how much longer the program will take.
        int nvox_remain = size_z * size_x * size_y;
        int counter = 0;
        int pref_ratio = 0;

        // --------------------------------------------------------------------
        // Smoothing loop
        cout << "  Smoothing with Gaussian..." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    counter++;
                    if ((counter * 100) / nvox_remain != pref_ratio) {
                        float c = (counter * 100) / nvox_remain;
                        cout << "\r  Progress: %" << c << flush;
                        pref_ratio = c;
                    }
                    for (int it = 0; it < size_t; ++it) {
                        int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                        *(nii_gaussw_data + voxel_i) = 0;

                        for (int it_i = max(0, it - vinc);
                             it_i < min(it + vinc + 1, size_t); ++it_i) {
                            int voxel_j = nxyz * it_i + nxy * iz + nx * iy + ix;

                            if (*(nii_data + voxel_j) != 0) {
                                float g = gaus(abs(it - it_i), gFWHM_val);
                                *(nii_smooth_data + voxel_i) +=
                                    *(nii_data + voxel_j) * g;
                                *(nii_gaussw_data + voxel_i) += g;
                            }
                        }
                    }
                }
            }
        }
        cout << endl;
    }

    // ========================================================================
    // Smoothing loop (for box-car smoothing)
    // ========================================================================
    if (do_box) {
        // float_kernel_size = 10;  // Corresponds to one voxel size.
        int vinc = bFWHM_val;  // Ignore if voxel is too far.
        cout << "  vinc " << vinc << endl;

        // Estimate how long this program will take.
        int nvox_remain = size_z * size_x * size_y;
        int counter = 0;
        int pref_ratio = 0;

        // --------------------------------------------------------------------
        // Smoothing loop //
        cout << "  Smoothing with box-car function." << endl;
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    counter++;
                    if ((counter * 100) / nvox_remain != pref_ratio) {
                        float c = (counter * 100) / nvox_remain;
                        cout << "\r  Progress: %" << c << flush;
                        pref_ratio = c;
                    }
                    for (int it = 0; it < size_t; ++it) {
                        int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                        *(nii_gaussw_data + voxel_i) = 0;

                        for (int it_i = max(0, it - vinc);
                             it_i < min(it + vinc + 1, size_t); ++it_i) {
                            int voxel_j = nxyz * it_i + nxy * iz + nx * iy + ix;

                            if (*(nii_data + voxel_j) != 0) {
                                *(nii_smooth_data + voxel_i) +=
                                    *(nii_data + voxel_j);
                                *(nii_gaussw_data + voxel_i) += 1;
                            }
                        }
                    }
                }
            }
        }
        cout << endl;
    }

    // ========================================================================
    // Correcting for edge error
    // ========================================================================
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_gaussw_data + i) != 0) {
            *(nii_smooth_data + i) /= *(nii_gaussw_data + i);
        }
    }

    nii_smooth->scl_slope = nii_input->scl_slope;
    if (nii_input->scl_inter != 0) {
        cout << "  ########################################## " << endl;
        cout << "  #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << "  ## The NIFTI scale factor is asymmetric ## " << endl;
        cout << "  ## Why would you do such a thing????    ## " << endl;
        cout << "  #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << "  ########################################## " << endl;
    }

    save_output_nifti(fin, "nii_smooth", nii_smooth, true);

    cout << "  Finished." << endl;
    return 0;
}
