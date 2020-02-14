
// Compile with with "make My_nii_read"
// Execute with ./My_nii_read -input input_example.nii -output output.nii -cutoff 3

#include <iostream>
#include <string>
#include <stdio.h>
#include "nifti2_io.h"
// #include "nifti2.h"
// #include "nifti1.h"
// #include "nifticdf.h"
// #include "nifti_tool.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

int show_help(void) {
    printf(
    "LN_FAsim: Short example of reading/writing NIfTI2.\n"
    "\n"
    "    This program is to demonstrate how to read a NIfTI-2 dataset.\n"
    "    Set output filenames and write a NIfTI-2 dataset, all via the\n"
    "    standard NIfTI C library.\n"
    "\n"
    "Usage:\n"
    "    LN_FAsim -input input_example.nii -output output.nii -cutoff 3 \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help.\n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -input  INFILE      : Specify input dataset.\n"
    "    -output OUTFILE     : Specify output dataset.\n"
    "    -verb LEVEL         : Set the verbose level to LEVEL.\n"
    "    -cutoff    value    : Set a cutof.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image *nim_input = NULL;
    char *fin = NULL, *fout = NULL;
    int cutoff, ac, disp_float_eg = 0;

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
            fin = argv[ac];  // no string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-verb")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -verb\n");
                return 2;
            }
            nifti_set_debug_level(atoi(argv[ac]));
        } else if (!strcmp(argv[ac], "-cutoff")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            cutoff = atoi(argv[ac]);  // no string copy, just pointer assignment
            cout << "Cutoff is " << cutoff << endl;
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

    cout << sizeSlice << " Slices | " <<  sizePhase << " PhaseSteps | " <<  sizeRead << " Read steps | " <<  nrep << " Timesteps "  << endl;

    if (!fout) {
        fprintf(stderr, "-- no output requested \n");
        return 0;
    }
    // Assign nifti_image fname/iname pair, based on output filename
    // (request to 'check' image and 'set_byte_order' here)
    if (nifti_set_filenames(nim_input, fout, 1, 1)) {
        return 1;
    }

    // Get access to data of nim_input
    float *nim_input_data = (float *) nim_input->data;

    // Allocating an additional nii
    nifti_image * nim_output1 = nifti_image_read(fin, 1);
    float *nim_output1_data = (float *) nim_output1->data;
    nim_output1->dim[4] = 1;
    // Changing according sizes nt etc.
    nifti_update_dims_from_array(nim_output1);

    // Allocating an additional nii
    nifti_image * nim_output2 = nifti_image_read(fin, 1);
    float *nim_output2_data = (float *) nim_output2->data;
    nim_output2->dim[4] = 1;
    // Changing according sizes nt etc.
    nifti_update_dims_from_array(nim_output2);

    // Get input as 4D array
    cout <<  "  Loading to array works 0 "  << endl;
    cout <<  "  Loading to array works 1 "  << endl;
    cout <<  "  Loading to array works 2"  << endl;

    // My individual analysis
    cout << "  Ok"  << endl;
    cout << "  Not ok"  << endl;

    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 1; ix < sizeRead - 1; ix = ix +2) {
                    *(nim_output1_data + nxyz * timestep + nxy * islice + nx * ix + iy) = -1 * *(nim_input_data + nxyz * timestep + nxy * islice + nx * ix + iy);
                }
            }
        }
    }
    for (int timestep = 0; timestep < nrep; ++timestep) {
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_output2_data + nxyz * timestep + nxy * islice + nx * ix + iy) = *(nim_input_data + nxyz * timestep + nxy * islice + nx * ix + iy) + 0.5* *(nim_input_data + nxyz * timestep + nxy * islice + nx * (ix - 1) + iy-1) + 0.5* *(nim_input_data + nxyz * timestep + nxy * islice + nx * (ix + 1) + iy);
                }
            }
        }
    }

    // Output file name
    const char *fout_3 = "zigzag.nii";
    if (nifti_set_filenames(nim_output1, fout_3 , 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_output1);

    // Output file name
    const char *fout_4 = "blurr.nii";
    if (nifti_set_filenames(nim_output2, fout_4 , 1, 1)) {
        return 1;
    }
    nifti_image_write(nim_output2);

    // Writing out input file
    // If we get here, write the output dataset

    // if (nifti_set_filenames(nim_input, fout, 1, 1)) {
    //     return 1;
    // }
    nifti_image_write(nim_input);
    // And clean up memory
    nifti_image_free(nim_output2);
    nifti_image_free(nim_output1);
    // nifti_image_free(nim_input);

    cout << "  Finished." << endl;
    return 0;
}
