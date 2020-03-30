#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_NOISE_KERNEL: \n"
    "         This program aims to estimate the functional point spread function\n"
    "         by looking at shared sources of variance across different directions.\n"
    "\n"
    "Usage:\n"
    "    LN_NOISE_KERNEL -input timeseries.nii  \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii) time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac;
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

    log_welcome("LN_NOISE_KERNEL");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);
    
    int kernal_size = 13; // This is the maximal number of layers. I don't know how to allocate it dynamically.
    float Nkernal[kernal_size][kernal_size][kernal_size] ; 
    int Number_of_averages = 0; 
    float Number_AVERAG[kernal_size][kernal_size][kernal_size] ; 

/////////////////////////////////
////////allokate and se zero ////
/////////////////////////////////

double vec_n[size_time] ;
double vec_nn[size_time] ;

    for(int timestep = 0; timestep < size_time ; timestep++) {
        vec_n[timestep] = 0; 
        vec_nn[timestep] = 0; 
    }


for(int i = 0; i < kernal_size; i++) {
    for(int j = 0; j < kernal_size ; j++) {
        for(int k = 0; k < kernal_size ; k++) {
            Nkernal[i][j][k] = 0; 
            Number_AVERAG[i][j][k] = 0; 
        }
    }
}


    // Allocate new nifti
    nifti_image* nii_kernel = nifti_copy_nim_info(nii);
    nii_kernel->nt = 1;
    nii_kernel->nx = kernal_size;
    nii_kernel->ny = kernal_size;
    nii_kernel->nz = kernal_size;
    nii_kernel->nvox =  kernal_size *  kernal_size *  kernal_size;
    nii_kernel->datatype = NIFTI_TYPE_FLOAT32;
    nii_kernel->nbyper = sizeof(float);
    nii_kernel->data = calloc(nii_kernel->nvox, nii_kernel->nbyper);
    float* nii_kernel_data = static_cast<float*>(nii_kernel->data);


    // ========================================================================
/* cout << "  Calculating skew, kurtosis, and autocorrelation..." << endl;

    double vec1[size_time];
    double vec2[size_time];

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it) {
                    vec1[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_skew_data + voxel_i) = ren_skew(vec1, size_time);
                *(nii_kurt_data + voxel_i) = ren_kurt(vec1, size_time);
                *(nii_autocorr_data + voxel_i) = ren_autocor(vec1, size_time);
            }
        }
    }

    save_output_nifti(fin, "skew", nii_skew, true);
    save_output_nifti(fin, "kurt", nii_kurt, true);
    save_output_nifti(fin, "autocorr", nii_autocorr, true);

    // ========================================================================
    cout << "  Calculating correlation with everything..." << endl;

    for (int it = 0; it < size_time; ++it) {
        vec1[it] = 0;
        vec2[it] = 0;
    }

    // Mean time course of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it) {
                    vec1[it] +=
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i)
                                            / nxyz);
                }
            }
        }
    }

    // Voxel-wise corelation to mean of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix <size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it)   {
                    vec2[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_conc_data + voxel_i) = ren_correl(vec1, vec2, size_time);
            }
        }
    }
    
*/ 

    nii_kernel->scl_slope = 1; // to make sure that the units are given in Pearson correlations (-1...1) 
    save_output_nifti(fin, "fPSF", nii_kernel, true);

    cout << "  Finished." << endl;
    return 0;
}
