


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_CORREL2FILES: Estimate the voxel wise correlation of two timeseries.\n"
    "\n"
    "    This program is motivated by Eli Merriam comparing in hunting down \n"
    "    voxels that out of phase for VASO and BOLD. \n"
    "\n"
    "Options:\n"
    "    -help  : Show this help.\n"
    "    -file1 : First time series.\n"
    "    -file2 : Second time series with should have the same dimensions \n"
    "             as first time series.\n"
    "    -output    : (Optional) Custom output name. \n"
    "                 including the path, if you want to write it as specific locations \n"
    "                 including the file extension: nii or nii.gz \n"
    "                 This will overwrite excisting files with the same name \n"
    "\n"
    "\n"
    "\n"
    "    an example application is mentioned on the blog post here: \n"
    "    http://layerfmri.com/QA \n"
    "\n"
    "Usage:\n"
    "    LN_CORREL2FILES -file1 file1.nii -file2 file2.nii \n"
    "                   \n"
    "    test application in the test_data folder would be:\n"
    "    ../LN_CORREL2FILES -file1 lo_Nulled_intemp.nii -file2 lo_BOLD_intemp.nii \n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ; 
    char *fin_1 = NULL, *fin_2 = NULL;
    int ac;
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-file1")) {
            cout << "Hello " << endl;
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file1\n");
                return 1;
            }
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-file2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file2\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-file1'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-file2'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
        return 2;
    }
    nifti_image*nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_CORREL2FILES");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    int size_z = nii1->nz;
    int size_x = nii1->nx;
    int size_y = nii1->ny;
    int size_time = nii1->nt;
    int nx = nii1->nx;
    int nxy = nii1->nx * nii1->ny;
    int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii1_temp = copy_nifti_as_float32(nii1);
    float* nii1_temp_data = static_cast<float*>(nii1_temp->data);
    nifti_image* nii2_temp = copy_nifti_as_float32(nii2);
    float* nii2_temp_data = static_cast<float*>(nii2_temp->data);

    // Allocate new nifti
    nifti_image *correl_file = nifti_copy_nim_info(nii1_temp);
    correl_file->nt = 1;
    correl_file->nvox = size_x * size_y * size_z;
    correl_file->data = calloc(correl_file->nvox, correl_file->nbyper);
    float *correl_file_data = static_cast<float*>(correl_file->data);
    // ========================================================================

    double vec1[size_time], vec2[size_time];
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it) {
                    int voxel_j = nxyz * it + nxy * iz + nx * iy + ix;
                    vec1[it] = *(nii1_temp_data + voxel_j);
                    vec2[it] = *(nii2_temp_data + voxel_j);
                }
                *(correl_file_data + voxel_i) =
                    static_cast<float>(ren_correl(vec1, vec2, size_time));
            }
        }
    }
    
    if (!use_outpath) fout = fin_1;
    save_output_nifti(fout, "correlated", correl_file, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
