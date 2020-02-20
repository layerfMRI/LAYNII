
#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_EXTREMETR: To find extreme the TR time point.\n"
    "\n"
    "    This program tries to find the maximal/minimal value of a time \n"
    "    series and writes out the timepoint of that TR vor every voxel.\n"
    "\n"
    "Usage:\n"
    "    LN_EXTREMETR -file file.nii \n"
    "\n"
    "Options:\n"
    "    -help : Show this help.\n"
    "    -file : Input time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    // nifti_image* nim_input=NULL;
    char *fin_1 = NULL;
    int ac;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file1\n");
                return 1;
            }
            fin_1 = argv[ac];  // Assign pointer, no string copy
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-file'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
      fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
      return 2;
    }

    log_welcome("LN_EXTREMETR");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    int size_x = nii1->nx;
    int size_y = nii1->ny;
    int size_z = nii1->nz;
    int size_t = nii1->nt;
    int nx = nii1->nx;
    int nxy = nii1->nx * nii1->ny;
    int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues

    nifti_image* nii1_temp = recreate_nii_with_float_datatype(nii1);
    float* nii1_temp_data = static_cast<float*>(nii1_temp->data);

    // ========================================================================

    nifti_image* max_file = nifti_copy_nim_info(nii1_temp);
    max_file->nt = 1;
    max_file->nvox = nii1_temp->nvox / size_t;
    max_file->datatype = NIFTI_TYPE_INT16;
    max_file->nbyper = sizeof(float);
    max_file->data = calloc(max_file->nvox, max_file->nbyper);
    short* max_file_data = (short*) max_file->data;

    nifti_image* min_file = nifti_copy_nim_info(nii1_temp);
    min_file->nt = 1;
    min_file->nvox = nii1_temp->nvox / size_t;
    min_file->datatype = NIFTI_TYPE_INT16;
    min_file->nbyper = sizeof(float);
    min_file->data = calloc(min_file->nvox, min_file->nbyper);
    short* min_file_data = (short*) min_file->data;

    double max_val = 0;
    double min_val = 1000000000;
    int TR_max = 0;
    int TR_min = 0;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                max_val = 0;
                min_val = 1000000000;
                for (int it = 0; it < size_t; ++it) {
                    if (*(nii1_temp_data + VOXEL_ID) > max_val) {
                        max_val = *(nii1_temp_data + VOXEL_ID);
                        TR_max = it;
                    }
                }
                for (int it = 0; it < size_t; ++it) {
                    if (*(nii1_temp_data + VOXEL_ID) < min_val) {
                        min_val = *(nii1_temp_data + VOXEL_ID);
                        TR_min = it;
                    }
                }
                *(min_file_data + VOXEL_ID_3D) = TR_min;
                *(max_file_data + VOXEL_ID_3D) = TR_max;
           }
        }
    }

    cout << "  Runing also until here 5... " << endl;

    string filename_1 = (string) (fin_1);
    string prefix_1 = "max_TR_";
    string outfilename_1 = prefix_1+filename_1;

    log_output(outfilename_1.c_str());
    const char *fout_1 = outfilename_1.c_str();
    if (nifti_set_filenames(max_file, fout_1 , 1, 1)) {
        return 1;
    }
    nifti_image_write(max_file);

    string prefix_2 = "min_TR_";
    string outfilename_2 = prefix_2+filename_1;
    log_output(outfilename_2.c_str());
    const char *fout_2 = outfilename_2.c_str();
    if (nifti_set_filenames(min_file, fout_2 , 1, 1)) {
        return 1;
    }
    nifti_image_write(min_file);

    // const char *fout_5 = "debug_ing.nii";
    // if (nifti_set_filenames(growfromWM0, fout_5 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(growfromWM0);

    // const char *fout_6 = "kootrGM.nii";
    // if (nifti_set_filenames(GMkoord2, fout_6 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(GMkoord2);

    // koord.autowrite("koordinaten.nii", wopts, &prot);
    cout << "  Finished." << endl;
    return 0;
}
