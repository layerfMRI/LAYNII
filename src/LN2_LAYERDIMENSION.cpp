#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_LAYERDIMENSION: This program switches the layer dimensions into nifti time dimension.\n"
    "                    This can be useful to browse layer profiles using the time course viewers\n"
    "                    of FSLEYES, AFNI, or miview. Furthermore, this is useful to execute\n"
    "                    time course analyses in the layer domain (e.g., ICA across layer profiles).\n"
    "\n"
    "Usage:\n"
    "    LN2_LAYERDIMENSION -values activation.nii -columns columns.nii -layers layers_equidist.nii -singleTR\n"
    "    ../LN2_LAYERDIMENSION -values lo_BOLD_act.nii -layers lo_layers.nii -columns lo_columns.nii \n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -values   : Nifti image with values that will be transformed into layer dimensions.\n"
    "                This is the contrast of interest, e.g. functional signal change.\n"
    "    -columns  : A 3D nifti file that contains columns as intager masks.\n"
    "                e.g. the output of LN2_COLUMNS or LN2_MULTILATERATE.\n"
    "    -layers   : A 3D nifti file that contains layers as intager masks.\n"
    "                For example LN2_LAYERS' output named 'layers'.\n"
    "                all three nii files above need to have the same spatial dimensions.\n"
    "    -output   : (Optional) Output filename, including .nii or\n"
    "                .nii.gz, and path if needed. Overwrites existing files.\n"
    "    -singleTR : flag to only look as the first time point of the value file.\n"
    "                default is ON.\n"
    "\n"
    "Notes:\n"
    "    - This does not refer to Dr. Strange's dimensions.\n"
    "           +-----------+ \n"
    "          /           /| \n"
    "         /           / | \n"
    "        /           /  | \n"
    "       +-----------+   | \n"
    "     L |           |   | \n"
    "     A |           |   + \n"
    "     Y |           |  / E\n"
    "     E |           | / M \n"
    "     R |           |/ I  \n"
    "       +-----------+ T   \n"
    "       SPACE             \n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL;
    int ac;
    bool mode_debug = false, mode_singleTR = true, use_outpath = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-columns")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -columns\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layers\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else if (!strcmp(argv[ac], "-singleTR")) {
            mode_singleTR = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-values'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-columns'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-layers'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }
    nii3 = nifti_image_read(fin3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }

    log_welcome("LN2_LAYERDIMENSION");
    log_nifti_descriptives(nii1); //values
    log_nifti_descriptives(nii2); //columns
    log_nifti_descriptives(nii3); //layers

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* layers = copy_nifti_as_int16(nii3);
    int16_t* layers_data = static_cast<int16_t*>(layers->data);
    nifti_image* columns = copy_nifti_as_int16(nii2);
    int16_t* columns_data = static_cast<int16_t*>(columns->data);

    // ========================================================================
    // Make sure there is nothing weird with the slope of the nii header
    // ========================================================================
    if(mode_debug){
       cout << "   Act  file has slope  " << nii1->scl_slope  << endl;
    }

    if (nii_input->scl_slope == 0) {
        cout << "   It seems like the slope of the value file is ZERO" << endl;
        cout << "   I am setting it to 1 instead " << endl;
        nii_input->scl_slope = 1;
    }

    // ========================================================================
    // Look how many layers and columns we have and allocate arrays accordingly
    // ========================================================================
    int nr_layers = 0;
    int nr_columns = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(layers_data + i) >= nr_layers){
            nr_layers = *(layers_data + i);
        }
        if (*(columns_data + i) >= nr_columns){
            nr_columns = *(columns_data + i);
        }
    }
    cout << "    There are " << nr_layers << " layers. " << endl;
    cout << "    There are " << nr_columns << " columns. " << endl << endl;

    // double numb_voxels[nr_layers][nr_columns];
    // double mean_val[nr_layers][nr_columns];
    std::vector<std::vector<double>> numb_voxels( nr_layers, std::vector<double>(nr_columns) );
    std::vector<std::vector<double>> mean_val( nr_layers,std::vector<double>(nr_columns) );

    for (int i = 0; i < nr_layers; i++) {
        for (int j = 0; j < nr_columns; j++) {
            mean_val   [i][j] = 0.;
            numb_voxels[i][j] = 0.;
        }
    }

    // ========================================================================
    // Prepare outputs
    // ========================================================================
    nifti_image* layerdim = copy_nifti_as_float32(nii_input);

    // Allocating new nifti for multi-dimensional images
    layerdim->datatype = NIFTI_TYPE_FLOAT32;
    layerdim->dim[0] = 4;  // For proper 4D nifti
    layerdim->dim[1] = nii_input->dim[1];
    layerdim->dim[2] = nii_input->dim[2];
    layerdim->dim[3] = nii_input->dim[3];
    layerdim->dim[4] = nr_layers;
    layerdim->nt = nr_layers;
    nifti_update_dims_from_array(layerdim);

    layerdim->nvox = nii_input->nvox * nr_layers;
    layerdim->nbyper = sizeof(float);
    layerdim->data = calloc(layerdim->nvox, layerdim->nbyper);
    layerdim->scl_slope = nii_input->scl_slope;
    layerdim->scl_inter = 0;
    float* layerdim_data = static_cast<float*>(layerdim->data);

    for (int voxi = 0; voxi < nr_voxels * nr_layers; voxi++) *(layerdim_data + voxi) = 0.0;

    // ========================================================================
    // Average within columns and layers
    // ========================================================================
    for (int i = 0; i != nr_voxels; ++i) {
        if ( *(columns_data + i) !=0 && *(layers_data + i) != 0 ){
            mean_val   [*(layers_data + i) -1 ][ *(columns_data + i) -1 ] += *(nii_input_data + i);
            numb_voxels[*(layers_data + i) -1 ][ *(columns_data + i) -1 ] += 1;
        }
    }

    for (int i = 0; i < nr_layers; i++) {
        for (int j = 0; j < nr_columns; j++) {
            if (numb_voxels[i][j] != 0){
                mean_val[i][j] /= (float)numb_voxels[i][j];
               // cout << "layer " << i+1 << " and column " << j+1 << " has value " << mean_val[i][j] << endl;
            }
        }
    }

    // ========================================================================
    // Fill average results into layer dimension file
    // ========================================================================
    for (int voxi = 0; voxi < nr_voxels; voxi++) {
        if ( *(columns_data + voxi) !=0 && *(layers_data + voxi) != 0){
            for (int l = 0; l < nr_layers; ++l) {
                *(layerdim_data + nr_voxels * l + voxi) = mean_val[ l ][ *(columns_data + voxi) - 1 ];
            }
            // *(layerdim_data + nr_voxels * 0 + voxi) = voxi;
        }
    }

    // ========================================================================
    // Write output
    // ========================================================================
    layerdim->pixdim[4] = 1/nr_layers; // in units of cortical depth.
    if (!use_outpath) fout = fin1;
    save_output_nifti(fout, "layerdim", layerdim, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
