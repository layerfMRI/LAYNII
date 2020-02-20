
// TODO(@Faruk): Compiles but gives error:
// `** ERROR: nifti_image_write_hdr_img: no image data`
// Also seems to give the same error in the old version. Check with Renzo.


#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_ZOOM: Reduces the nii file to only include the rectangular space \n"
    "         that is different from zero.\n"
    "\n"
    "Usage: \n"
    "    LN_ZOOM -mask mask.nii -input file_to_be_zoomed.nii \n"
    "\n"
    "Options:\n"
    "\n"
    "    -help  : Show this help.\n"
    "    -input : Nifti (.nii) file that should be zoomed (e.g. with \n"
    "             multiple time points).\n"
    "    -mask  : Nifti (.nii) file that determines the region of interest\n"
    "             (e.g. the layer mask with one time point).\n"
    "\n"
    "Note: Written for Insub. \n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    // nifti_image* nim_input=NULL;
    char* fin_1 = NULL, *fin_2 = NULL;
    int ac, disp_float_eg = 0, shift = 0;
    int trialdur = 0;
    if (argc < 2) return show_help();   // Typing '-help' is sooo much work

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin_1 = argv[ac];  // no string copy, just pointer assignment
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            fin_2 = argv[ac];  // no string copy, just pointer assignment
        } else {
            fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-mask'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_ZOOM");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);  // mask

    // Get dimensions of input
    const int size_z = nii1->nz;
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_t = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues

    nifti_image* nim_file_1 = recreate_nii_with_float_datatype(nii1);
    float* nii1_data = static_cast<float*>(nim_file_1->data);

    nifti_image* nim_file_2 = recreate_nii_with_float_datatype(nii2);
    float* nii2_data = static_cast<float*>(nim_file_2->data);

    // ========================================================================

    // Initialize min-max inversely to increas-decrease in the loop
    int min_z = 10000, max_z = 0;
    int min_x = 10000, max_x = 0;
    int min_y = 10000, max_y = 0;
    float mask_val = 0;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                mask_val = *(nii2_data + VOXEL_ID_3D);
                if (mask_val > 0 && ix > max_x) max_x = ix;
                if (mask_val > 0 && iy > max_y) max_y = iy;
                if (mask_val > 0 && iz > max_z) max_z = iz;
                if (mask_val > 0 && ix < min_x) min_x = ix;
                if (mask_val > 0 && iy < min_y) min_y = iy;
                if (mask_val > 0 && iz < min_z) min_z = iz;
            }
        }
    }
    cout << "  x range is " << min_x << "-" << max_x << endl;
    cout << "  y range is " << min_y << "-" << max_y << endl;
    cout << "  z range is " << min_z << "-" << max_z << endl;
    int zoomed_z_size = max_z - min_z + 1;
    int zoomed_x_size = max_x - min_x + 1;
    int zoomed_y_size = max_y - min_y + 1;

    // Handle zoomed nifti
    nifti_image* zoomed_file = nifti_copy_nim_info(nim_file_1);
    // zoomed_file->nt = 7;
    zoomed_file->nz = zoomed_z_size;
    zoomed_file->nx = zoomed_x_size;
    zoomed_file->ny = zoomed_y_size;

    zoomed_file->nvox = (nim_file_1->nvox * zoomed_x_size * zoomed_y_size * zoomed_z_size) / (nii1->nx * nii1->ny * nii1->nz);
    zoomed_file->datatype = NIFTI_TYPE_FLOAT32;
    zoomed_file->nbyper = sizeof(float);
    zoomed_file->data = calloc(zoomed_file->nvox, zoomed_file->nbyper);
    float* zoomed_file_data = static_cast<float*>(zoomed_file->data);

    cout << "  Reduction " << (nim_file_1->nvox * zoomed_x_size * zoomed_y_size * zoomed_z_size) / (nii1->nx * nii1->ny * nii1->nz) << endl;

    const int nx_2 = zoomed_file->nx;
    const int nxy_2 = zoomed_file->nx * zoomed_file->ny;
    const int nxyz_2 = zoomed_file->nx * zoomed_file->ny * zoomed_file->nz;

    for (int it = 0; it < size_t; ++it) {
        for (int iz = min_z; iz < max_z+1; ++iz) {
            for (int iy = min_y; iy < max_y+1; ++iy) {
                for (int ix = min_x; ix < max_x+1; ++ix) {
                    *(zoomed_file_data + nxyz_2 * it + nxy_2 * (iz - min_z) + nx_2 * (iy - min_y) + (ix - min_x))
                        = *(nii1_data + nxyz * it + nxy * iz + nx * iy + ix);
                }
            }
        }
    }

    string prefix = "zoomed_";
    string filename_1 = (string) (fin_1);
    string outfilename = prefix+filename_1;
    const char* fout_1 = outfilename.c_str();
    if (nifti_set_filenames(zoomed_file, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(zoomed_file);
    log_output(outfilename.c_str());

    return 0;
}
