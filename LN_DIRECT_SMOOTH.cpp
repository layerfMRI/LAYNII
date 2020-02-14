
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

using namespace std;

#define PI 3.14159265;

// #include "utils.hpp"

int show_help(void) {
    printf(
    "LN_DIRECT_SMOOTH : Smoothing in specific directions only. This program \n"
    "                   smooths data within layer or columns. In order to \n"
    "                   avoid smoothing across masks a crawler smooths only \n"
    "                   across connected voxels.\n"
    "\n"
    "Usage:\n"
    "    LN_DIRECT_SMOOTH -input activity_map.nii -FWHM 1 -direction 1 \n"
    "    LN_DIRECT_SMOOTH -input bico_VASO.Mean.nii -FWHM 0.5 -direction 3 -laurenzian \n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti (.nii) file that will be smoothed. It \n"
    "                     should have same dimensions as layer file.\n"
    "    -FWHM             : Amount of smoothing in units of voxels.\n"
    "    -laurenzian    : Use Laurenzian smoothing. Default is Gaussian \n"
    "                   : only for division images.\n"
    "    -direction     : Axis of smoothing. 1 for x, 2 for y or 3 for z. \n"
    "    -Anonymous_sri : You know what you did (no FWHM).\n"
    "\n"
    "Note: This program supports INT16, INT32 and FLOAT32.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * fout=NULL, * finfi=NULL;
    int ac, direction_i=0, do_laurenz = 0, do_gauss = 0, do_sri = 0;
    float FWHM_val=10, strength = 1;
    if (argc < 3) {
        return show_help();  // Typing '-help' is sooo much work
    }
    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
            fprintf(stderr, "** missing argument for -FWHM\n");
            return 1;
            }
            FWHM_val = atof(argv[ac]);  // Pointer assignment, no string copy.
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            finfi = argv[ac];  // Pointer assignment, no string copy.
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            direction_i = atoi(argv[ac]);  // Pointer assignment, no string copy.
        } else if (!strcmp(argv[ac], "-laurenzian")) {
            do_laurenz = 1;
            fprintf(stderr, "Use Laurenzian smoothing instead of Gaussian.");
        } else if (!strcmp(argv[ac], "-Anonymous_sri")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Anonymous_sri\n");
            // return 1;
            }
            do_sri = 1;
            strength = atof(argv[ac]);  // Pointer assignment, no string copy.
            fprintf(stderr, "Yes Sri I am doing you ");
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finfi) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
    if (!nim_inputfi) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", finfi);
        return 2;
    }
    if (direction_i == 0) {
        fprintf(stderr, "** failed to read direction '%i'\n", direction_i);
        return 2;
    }
    // Get dimensions of input
    int sizeSlice = nim_inputfi->nz;
    int sizePhase = nim_inputfi->nx;
    int sizeRead = nim_inputfi->ny;
    int nrep =  nim_inputfi->nt;
    int nx =  nim_inputfi->nx;
    int nxy = nim_inputfi->nx * nim_inputfi->ny;
    int nxyz = nim_inputfi->nx * nim_inputfi->ny * nim_inputfi->nz;
    float dX =  1;  // nim_inputfi->pixdim[1];
    float dY =  1;  // nim_inputfi->pixdim[2];
    float dZ =  1;  // nim_inputfi->pixdim[3];

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_inputf = nifti_copy_nim_info(nim_inputfi);
    nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
    nim_inputf->nbyper = sizeof(float);
    nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
    float *nim_inputf_data = (float *) nim_inputf->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////

    if (nim_inputfi->datatype == NIFTI_TYPE_FLOAT32 ||  nim_inputfi->datatype == NIFTI_TYPE_INT32) {
        float *nim_inputfi_data = (float *) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz *it + nxy*islice + nx*ix + iy) = (float) (*(nim_inputfi_data + nxyz *it + nxy*islice + nx*ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == NIFTI_TYPE_INT16 || nim_inputfi->datatype == DT_UINT16) {
        short *nim_inputfi_data = (short *) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix=0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz *it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == DT_FLOAT64 || nim_inputfi->datatype == NIFTI_TYPE_FLOAT64) {
        double *nim_inputfi_data = (double *) nim_inputfi->data;
        for (int it = 0; it< nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    cout << sizeSlice << " Slices    " << sizePhase << " PhaseSteps     " << sizeRead << " Read steps    " << nrep << " timesteps " << endl;
    cout << "  Voxel size = " << dX << " x " << dY << " x " << dZ << endl;
    cout << "  Datatype 1 = " << nim_inputf->datatype << endl;

    /////////////////////////////////////
    // Make allocating necessary files //
    /////////////////////////////////////
    nifti_image * smoothed = nifti_copy_nim_info(nim_inputf);
    nifti_image * gausweight = nifti_copy_nim_info(nim_inputf);
    // nifti_image * layer = nifti_copy_nim_info(nim_input);
    // nifti_image * leak_layer = nifti_copy_nim_info(nim_input);

    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    // layer->datatype = NIFTI_TYPE_FLOAT32;
    // leak_layer->datatype = NIFTI_TYPE_FLOAT32;

    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    // layer->nbyper = sizeof(float);
    // leak_layer->nbyper = sizeof(float);

    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
    // layer->data = calloc(layer->nvox, layer->nbyper);
    // leak_layer->data = calloc(leak_layer->nvox, leak_layer->nbyper);

    float *smoothed_data = (float *) smoothed->data;
    float *gausweight_data = (float *) gausweight->data;
    // float *layer_data = (float *) layer->data;
    // float *leak_layer_data = (float *) leak_layer->data;
    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float gaus(float distance, float sigma);
    float laur(float distance, float sigma);
    float ASLFt(float distance, float strength);

    cout << "  Debug 2 " << endl;

    // Float kernal_size = 10;  // corresponds to one voxel sice.
    int vinc = max(1., 2. * FWHM_val / dX);  // Ignore far voxels
    float dist_i = 0.;
    cout << "  vinc " << vinc<< endl;
    cout << "  FWHM_val " << FWHM_val<< endl;

    if (do_sri == 0 && do_laurenz == 0) {
        do_gauss = 1;
    }
    if (do_sri == 1) {
        cout << "  No FWHM_val is used. Strength is " << strength<< endl;
    }
    // cout << " No FWHM_val is used. Strength is " << strength<< endl;

    ////////////////////
    // Smoothing loop //
    ////////////////////
    // cout << " DEBUG " << dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl;
    cout << "  Smoothing in dimension " << direction_i << endl;
    cout << "  Timestep... " << flush;

    for (int i = 0; i < 5;  i++) {
        cout << "  " << ASLFt(i, strength) << endl;
    }

    if (direction_i == 1) {
        for (int it = 0; it < nrep; ++it) {
            cout << "  " << it << "..." << flush;
            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead - 0; ++ix) {
                        *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                        //*(smoothed_data + nxy * iz + nx * ix + iy) = 0;

                        for (int ix_i = max(0, ix - vinc); ix_i < min(ix + vinc + 1, sizeRead); ++ix_i) {
                            if (*(nim_inputf_data + nxyz * it + nxy * iz + nx * ix_i + iy) != 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy, (float)iz, dX, dY, dZ);

                                if (do_gauss == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix_i + iy) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz *it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix_i + iy) * laur(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz *it + nxy*iz + nx*ix + iy) + laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix_i + iy) * ASLFt(dist_i, strength);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + ASLFt(dist_i, strength);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (direction_i == 2) {
        for (int it = 0; it < nrep; ++it) {
            cout << "  " << it << "..." << flush;
            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead - 0; ++ix) {
                        *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                        //*(smoothed_data + nxy * iz + nx * ix + iy) = 0;

                        for (int iy_i = max(0, iy - vinc); iy_i < min(iy + vinc + 1, sizePhase); ++iy_i) {
                            if (*(nim_inputf_data + nxyz * it + nxy * iz + nx * ix + iy_i) != 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy_i, (float)iz, dX, dY, dZ);

                                if (do_gauss == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix + iy_i) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix + iy_i) * laur(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz + nx * ix + iy_i) * ASLFt(dist_i, strength);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + ASLFt(dist_i, strength);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (direction_i == 3) {
        for (int it = 0; it < nrep; ++it) {
            cout << "  " << it << "..." << flush;
            for (int iz = 0; iz < sizeSlice; ++iz) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead-0; ++ix) {
                        *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                        // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                        for (int iz_i=max(0,iz-vinc); iz_i<min(iz+vinc+1,sizeSlice); ++iz_i) {
                            if (*(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) != 0) {

                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy, (float)iz_i, dX, dY, dZ);

                                if (*(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) == 0) {
                                    cout << *(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) << endl;
                                }

                                if (do_gauss == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) * gaus(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz *it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) * laur(dist_i, FWHM_val);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix + iy) * ASLFt(dist_i, strength);
                                    *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + ASLFt(dist_i, strength);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << endl;

    ///////////////////////////////
    // Correcting for edge error //
    ///////////////////////////////
    for (int it = 0; it < nrep; ++it) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) / *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    cout << "  Running also until here  " << endl;

    smoothed->scl_slope =  nim_inputfi->scl_slope;

    if (nim_inputfi->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    // output file name
    // const char *fout_4 ="leaky_layers.nii";
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) return 1;
    // nifti_image_write(leak_layer);
    //
    // const char *fout_5 = "input_file.nii";
    // if (nifti_set_filenames(nim_inputf, fout_5, 1, 1)) return 1;
    // nifti_image_write(nim_inputf);
    //
    // const char *fout_2 = "mask.nii";
    // if (nifti_set_filenames(nim_mask, fout_2, 1, 1)) return 1;
    // nifti_image_write(nim_mask);

    string prefix = "smoothed_";
    string filename = (string) (finfi);
    string outfilename = prefix+filename;

    cout << "  Writing as = " << outfilename.c_str() << endl;  // finfi is: char *

    const char *fout_1 = outfilename.c_str();
    if (nifti_set_filenames(smoothed, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(smoothed);

    // const char *fout_1 = "layer.nii";
    // if (nifti_set_filenames(layer, fout_1, 1, 1)) return 1;
    // nifti_image_write(layer);

    cout << "  Finished." << endl;
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt((x1 - x2) * (x1-x2) * dX * dX
                + (y1 - y2) * (y1 - y2) * dY * dY
                + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592)) * exp (-0.5 * distance * distance / (sigma * sigma));
}

float laur (float distance, float sigma) {
    // Note: For consistency with Gaus's sigma, I am using a scaled version of
    // the FWHM.
    // sigma = sigma / sqrt(2 * log (2));
    float result = 0;
    if ((int)distance%2 == 0) {
        //sigma = 2 * sigma;
        result = 1 * 1. / (3.141592 * sigma) * 1 /(1 + distance * distance / (sigma * sigma));
    }
    if ((int)distance%2 == 1) {
        //sigma = 2 * sigma;
        result = 1.6 * 1. / (3.141592 * sigma) * 1 / (1 + distance * distance / (sigma * sigma));
    }
return result;
// return 1 / 3.141592 * 1/ ((distance) * (distance) + (0.5 * sigma) * (0.5 * sigma));
}

float ASLFt(float distance, float strength) {
    float value;
    if (distance >= 4.0) value = 0.0;
    if (distance <= 4.0) value = 0.0245;
    if (distance <= 3.0) value = 0.024;
    if (distance <= 2.0) value = 0.06;
    if (distance <= 1.0) value = 0.18;
    if (distance == 0.0) value = 1.0;
    return value * strength;
}
