#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_CORREL2FILES: Estimate the voxel wise correlation of two timeseries.\n"
    "\n"
    "Usage:\n"
    "    LN_CORREL2FILES -file1 file1.nii -file2 file2.nii \n"
    "    ../LN_CORREL2FILES -file1 lo_Nulled_intemp.nii -file2 lo_BOLD_intemp.nii \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -file1  : First time series.\n"
    "    -file2  : Second time series with should have the same dimensions \n"
    "              as first time series.\n"
    "    -output : (Optional) Output filename, including .nii or\n"
    "              .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n"
    "Notes:\n"
    "    - This program is used for hunting voxels that are out of phase in\n"
    "      nulled and BOLD time series.\n"
    "    - An example application is mentioned in the blog post here:\n"
    "      <http://layerfmri.com/QA> \n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
    char *fout = NULL, *fin_1 = NULL, *fin_2 = NULL;
    int ac;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-file1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file1\n");
                return 1;
            }
            fin_1 = argv[ac];
            fout = argv[ac];
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
    int64_t size_z = nii1->nz;
    int64_t size_x = nii1->nx;
    int64_t size_y = nii1->ny;
    int64_t size_time = nii1->nt;
    int64_t nr_voxels = size_z * size_y * size_x;

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
    std::vector<double> vec1(size_time);
    std::vector<double> vec2(size_time);
    int64_t ix, iy, iz;
    for (int64_t i = 0; i < nr_voxels; ++i) {
        tie(ix, iy, iz) = ind2sub_3D_64(i, size_x, size_y);
        for (int64_t it = 0; it < size_time; ++it) {
            int64_t j = nr_voxels*it + size_x*size_y*iz + size_x*iy + ix;
            vec1[it] = *(nii1_temp_data + j);
            vec2[it] = *(nii2_temp_data + j);
        }
        *(correl_file_data + i) = static_cast<float>(ren_correl(vec1.data(), vec2.data(), size_time));
    }

    save_output_nifti(fout, "correlated", correl_file, true);

    cout << "  Finished." << endl;
    return 0;
}
