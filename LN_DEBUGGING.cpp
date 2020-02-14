
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
    "LN_DEBUGGINH: Short example of Layering.\n"
    "\n"
    "    This program demonstrates how to read a NIfTI-2 dataset.\n"
    "    Set output filenames and write a NIfTI-2 dataset, all via the\n"
    "    standard NIfTI C library.\n"
    "\n"
    "Usage:\n"
    "    LN_DEBUGGING -rim rim.nii \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help.\n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -rim border         : Specify input dataset.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image *nim_input = NULL;
    char *fin = NULL, *fout = NULL;
    int ac, disp_float_eg = 0;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // no string copy, just pointer assignment
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }
    // Read input dataset, including data
    nim_input = nifti_image_read(fin, 1);

    // Get dimensions of input
    int sizeSlice = nim_input->nz;
    int sizePhase = nim_input->nx;
    int sizeRead = nim_input->ny;
    int nrep = nim_input->nt;
    int nx = nim_input->nx;
    int nxy = nim_input->nx * nim_input->ny;
    int nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    cout << sizeSlice << " Slices | " << sizePhase << " PhaseSteps | " << sizeRead << " Read steps | " << nrep << " Timesteps " << endl;

    // Get access to data of nim_input
    nim_input = nifti_image_read(fin, 1);
    short *nim_input_data = (short *) nim_input->data;

    nifti_image * nim_2 = nifti_image_read(fin, 1);
    cout << "  1st datatype 1 " << nim_input->datatype << endl;
    // short  *nim_2_data = (short *) nim_2->data;

    nifti_image * nim_5 = nifti_image_read(fin, 1);
    nim_5->datatype = NIFTI_TYPE_FLOAT32;
    nim_5->nbyper = sizeof(float);
    nim_5->data = calloc(nim_5->nvox, nim_5->nbyper);
    float *nim_5_data = (float *) nim_5->data;

    if (nim_input->datatype == NIFTI_TYPE_FLOAT32) {
        float *nim_2_data = (float *) nim_2->data;
        cout << "  Datatype 1 " << nim_2->datatype << endl;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_5_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_2_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_input->datatype == NIFTI_TYPE_FLOAT32) {
        float  *nim_2_data = (float *) nim_2->data;
        cout << "  Datatype 1 " << nim_2->datatype << endl;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_5_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_2_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_input->datatype == NIFTI_TYPE_INT16) {
        short  *nim_2_data = (short *) nim_2->data;
        cout << "  Datatype 1 " << nim_2->datatype << endl;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_5_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_2_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_input->datatype == NIFTI_TYPE_INT32) {
        int *nim_2_data = (int *) nim_2->data;
        cout << "  Datatype 1 " << nim_2->datatype << endl;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_5_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_2_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    // Rick comments
    // nifti_copy_nim_info(). It will return with data == NULL.
    // If you need the data allocated, memory use would not change once you do so.
    // There is also nifti_make_new_nim()
    // nifti_image * nifti_make_new_nim(const int64_t dims[], int datatype, int data_fill)

    // for (int timestep = 0; timestep < 100; ++timestep) {
    //     cout << timestep << "  Timestep " << endl;
    //     for (int islice = 0; islice < sizeSlice; ++islice) {
    //         for (int iy = 0; iy < sizePhase; ++iy) {
    //             for (int ix = 0; ix < sizeRead; ++ix) {
    //                 // cout << islice << " Slices " << iy << " | PhaseSteps " << ix << " | Read steps " << endl;
    //               *(nim_3_data + nxyz * timestep + nxy * islice + nx * ix + iy) = 1;
    //             }
    //         }
    //     }
    // }

    cout << "  Runing also until here 5... " << endl;

     // Output file name
    const char  *fout_5 = "nim_5.nii";
    if (nifti_set_filenames(nim_5, fout_5 , 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_5);

    const char  *fout_2 = "nim_2.nii";
    if (nifti_set_filenames(nim_2, fout_2 , 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_2);

    cout << "  Finished." << endl;
    return 0;
  }
