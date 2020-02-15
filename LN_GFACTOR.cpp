
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "nifti2_io.h"
// #include "nifti2.h"
// #include "nifti1.h"
// #include "nifticdf.h"
// #include "nifti_tool.h"
// #include < gsl/gsl_multifit.h>
// #include < gsl/gsl_statistics_double.h>
using namespace std;

#define PI 3.14159265;

// #include "utils.hpp"

int N_rand = 1000;
double_t lower = -10;
double_t upper = 10;

double verteilung(double x);
typedef double ( * Functions)(double);
Functions pFunc = verteilung;

double_t arb_pdf_num(int N_rand, double ( * pFunc)(double), double_t lower,
                     double_t upper);
double adjusted_rand_numbers(double mean, double stdev, double value);

int show_help(void) {
    printf(
        "LN_GFACTOR: Simulating where the g-factor penalty would be largest.\n"
        "\n"
        "Usage:\n"
        "    LN_GFACTOR -input MEAN.nii -output GFACTOR.nii -variance 1 -direction 1 -grappa 2 -cutoff 150 \n"
        "\n"
        "Options:\n"
        "    -help               : Show this help.\n"
        "    -disp_float_example : Show some voxel's data.\n"
        "    -input  INFILE      : Specify input dataset.\n"
        "    -output OUTFILE     : Specify output dataset.\n"
        "    -variance value     : How much noise there will be.\n"
        "    -direction          : Phase encoding direction [0=x, 1=y, 2=z].\n"
        "    -grappa             : TODO(Renzo): Missing description."
        "    -cutoff value       : Value to seperate noise from signal.\n"
        "    -verb LEVEL         : Set the verbose level to LEVEL.\n"
        "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image * nim_input = NULL;
    char * fin = NULL, * fout = NULL;
    int grappa_int, direction_int, ac, disp_float_eg = 0;
    float cutoff, variance_val;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-disp_float_example")) {
            disp_float_eg = 1;
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-variance")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            variance_val = atof(argv[ac]);  // Assign pointer, no string copy
            cout << "Variance  = " << variance_val << endl;
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            direction_int = atof(argv[ac]);  // Assign pointer, no string copy
            cout << "Direction = " << direction_int << endl;
        } else if (!strcmp(argv[ac], "-grappa")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -grappa\n");
                return 1;
            }
            grappa_int = atof(argv[ac]);  // Assign pointer, no string copy
            cout << "GRAPPA    = " << grappa_int << endl;
        } else if (!strcmp(argv[ac], "-cutoff")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -cutoff\n");
                return 1;
            }
            cutoff = atof(argv[ac]);  // Assign pointer, no string copy
            cout << "Cut-off   = " << cutoff << endl;
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
    nim_input = nifti_image_read(fin, 1);
    if (!nim_input) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    // Get dimensions of input
    int sizeSlice = nim_input->nz;
    int sizePhase = nim_input->nx;
    int sizeRead = nim_input->ny;
    int nrep = nim_input->nt;
    int nx = nim_input->nx;
    int nxy = nim_input->nx * nim_input->ny;
    int nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    cout << sizeSlice << " Slices | " << sizePhase << " Phase_steps | " << sizeRead << " Read_steps | " << nrep << " Time_steps " << endl;
    // cout << sizePhase %grappa_int << " mod    " << endl;

    if (direction_int == 0) {
        sizePhase = sizePhase - (sizePhase %grappa_int);
    }
    if (direction_int == 1) {
        sizeRead = sizeRead  - (sizeRead  %grappa_int);
    }
    if (direction_int == 2) {
        sizeSlice = sizeSlice - (sizeSlice %grappa_int);
    }

    if (!fout) {
        fprintf(stderr, "-- no output requested. \n");
        return 0;
    }

    // Assign nifti_image fname/iname pair, based on output filename
    // (request to 'check' image and 'set_byte_order' here)
    if (nifti_set_filenames(nim_input, fout, 1, 1)) {
        return 1;
    }

    // Get access to data of nim_input
    float * nim_input_data = (float * ) nim_input->data;

    // Allocating additional images
    nifti_image * gfactormap = nifti_copy_nim_info(nim_input);
    gfactormap->datatype = NIFTI_TYPE_FLOAT32;
    gfactormap->nbyper = sizeof(float);
    gfactormap->data = calloc(gfactormap->nvox, gfactormap->nbyper);
    float * gfactormap_data = (float * ) gfactormap->data;

    // allocating additional images
    nifti_image * binary = nifti_copy_nim_info(nim_input);
    binary->datatype = NIFTI_TYPE_FLOAT32;
    binary->nbyper = sizeof(float);
    binary->data = calloc(binary->nvox, binary->nbyper);
    float * binary_data = (float * ) binary->data;

    // noise_image
    nifti_image * noise_image = nifti_copy_nim_info(nim_input);
    noise_image->datatype = NIFTI_TYPE_FLOAT32;
    noise_image->nbyper = sizeof(float);
    noise_image->data = calloc(noise_image->nvox, noise_image->nbyper);
    float * noise_image_data = (float * ) noise_image->data;

    // for (int timestep = 0; timestep < nrep; ++timestep) {
    //     for (int islice = 0; islice < sizeSlice; ++islice) {
    //         for (int iy = 0; iy < sizePhase; ++iy) {
    //             for (int ix = 0; ix < sizeRead; ++ix) {
    //                 *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(nim_input_data + nxyz * timestep + nxy * islice + nx * ix + iy) + adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper));
    //                 //cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
    //             }
    //         }
    //     }
    // }

    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (cutoff <= *(nim_input_data + nxyz * timestep + nxy * islice + nx * ix + iy)) {
                        *(binary_data + nxyz * timestep + nxy * islice + nx * ix + iy) = 1;
                    } else {
                        *(binary_data + nxyz * timestep + nxy * islice + nx * ix + iy) = 0;
                    }
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }

    const char * fout_7 = "binary.nii";
    if (nifti_set_filenames(binary, fout_7, 1, 1)) {
        return 1;
    }
    nifti_image_write(binary);

    // Setting initial condition
    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = 0.;
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }

    ////////////////////////////////////
    // Adding additional GRAPPA noise //
    ////////////////////////////////////

    if (direction_int == 0) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase/grappa_int; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) + *(binary_data + nxyz * timestep + nxy * islice + nx * ix + iy+ grappa_seg * sizePhase/grappa_int)/grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg=1; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase/grappa_int; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy+ grappa_seg * sizePhase / grappa_int) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    if (direction_int == 1) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead/grappa_int; ++ix) {
                            *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) + *(binary_data + nxyz * timestep + nxy * islice + nx * (ix+ grappa_seg * sizeRead/grappa_int) + iy)/grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead/grappa_int; ++ix) {
                            *(gfactormap_data + nxyz * timestep + nxy * islice + nx * (ix+ grappa_seg * sizeRead/grappa_int) + iy) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    if (direction_int == 2) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice/grappa_int; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) + *(binary_data + nxyz * timestep + nxy * (islice+ grappa_seg * sizeSlice/grappa_int) + nx * ix + iy)/grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int timestep = 0; timestep < nrep; ++timestep) {
                for (int islice = 0; islice < sizeSlice/grappa_int; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(gfactormap_data + nxyz * timestep +nxy * (islice+ grappa_seg * sizeSlice/grappa_int) + nx * ix + iy) = *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(binary_data + nxyz * timestep + nxy * islice + nx * ix + iy) **(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy);
                }
            }
        }
    }
    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(noise_image_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(nim_input_data + nxyz * timestep + nxy * islice + nx * ix + iy) + *(gfactormap_data + nxyz * timestep + nxy * islice + nx * ix + iy) * cutoff * adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper));
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }

    // Output file name
    if (nifti_set_filenames(gfactormap, fout, 1, 1)) {
        return 1;
    }
    nifti_image_write(gfactormap);

    const char * fout_2 = "amplified_GRAPPA.nii";
    if (nifti_set_filenames(noise_image, fout_2, 1, 1)) {
        return 1;
    }
    nifti_image_write(noise_image);

    // Writing out input file
    // if we get here, write the output dataset

    // if (nifti_set_filenames(nim_input, fout , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(nim_input);
    //  and clean up memory
    nifti_image_free(gfactormap);
    // nifti_image_free(nim_input);
    cout << "  Finished." << endl;
    return 0;
}

// Gauss lower = -5, upper = 5
double verteilung(double z) {
    return exp(-z * z / (2.)) * 1. / sqrt(2. * M_PI);
}

double_t arb_pdf_num(int N_rand, double ( * pFunc)(double), double_t lower,
                     double_t upper) {
    double_t binwidth = (upper - lower)/(double_t)N_rand;
    double_t integral = 0.0;
    double_t rand_num = rand() / (double_t)RAND_MAX;
    int i;

    for (i = 0; integral < rand_num; i++) {
        integral += pFunc(lower + (double_t) i * binwidth) * binwidth;
        if ((lower + (double_t) i * binwidth) > upper) {
            cout << "  Upper limit, vielleicht sollte da limit angepasst werden." << i << endl;
            return lower + (double_t) i * binwidth;
        }
    }
    return lower + (double_t) i * binwidth;
}

double adjusted_rand_numbers(double mean, double stdev, double value) {
    return value * stdev + mean;
}
