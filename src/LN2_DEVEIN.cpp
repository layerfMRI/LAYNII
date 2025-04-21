#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN2_DEVEIN : Mitigates venous signal in one of three different ways.\n"
    "\n"
    "    1. Estimates microvascular component in layer-fMRI GE-BOLD. It does so by using an estimate\n"
    "       of the local macrovascular blood volume amplitude of lod frequencies (ALF).\n"
    "    2. Linearly scales back input values using layers information.\n"
    "    3. Scales input values following CBV.\n"
    "\n"
    "Usage:\n"
    "    LN2_DEVEIN -input lo_BOLD_act.nii -layer_file lo_layers.nii -linear\n"
    "    LN2_DEVEIN -input lo_BOLD_act.nii -layer_file lo_layers.nii -column_file lo_columns.nii -ALF lo_ALF.nii\n"
    "    ../LN2_DEVEIN -layer_file lo_layers.nii -column_file lo_columns.nii -input lo_BOLD_act.nii -ALF lo_ALF.nii\n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -input       : BOLD file that should be corrected from macrovascular contaminations. This means\n"
    "                   that there will be no deconvolution, just scaling of the baseline venous CBV.\n"
    "                   This can be a time series or an activity map (not z-scores though). \n"
    "    -layer_file  : Nifti (.nii) file that contains layers. Layers close to white matter\n"
    "                   should be denoted by smaller numbers.\n"
    "    -column_file : Nifti (.nii) file that contains columns.\n"
    "    -ALF         : (Optional) File with estimates of amplitude of low frequencies (ALF)\n"
    "                   as a correlate to venous CBV. \n"
    "    -CBV         : (Optional) For CBV scaling \n"
    "    -linear      : (Optional) For linear scaling, this option should selected used, when no \n"
    "                   ALF file is given. \n"
    "    -lambda      : (Optional) For peak to tail ratio. Default is 0.25\n"
    "                   from Markuerkiaga et al. 2016, Fig. 5B, at 7T.\n"
    "    -output      : (Optional) Output filename, including .nii or\n"
    "                   .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n"
    "Notes:\n"
    "    - [On lambda parameter]: If you assume your cerebral blood flow (CBF) is exceptionally low,\n"
    "       use 0.3. If you assume your CBV is exceptionally high, use 0.2. If you use 3T instead of\n"
    "       7T, use 20 percent larger values. Larger values will result in stronger deconvolution.\n"
    "       These variations will not affect the resulting profiles too much. Therefore only change\n"
    "       this parameter if it is absolutely needed.\n"
    "    - This program is described in more depth in this blog post:\n"
    "       <https://layerfmri.com/devein> \n"
    "    - The leakage model is described in:\n"
    "       Markuerkiaga et al. (2016) <<https://doi.org/10.1016/j.neuroimage.2016.02.073>>\n"
    "       Havlicek and Uludag (2019) <<https://doi.org/10.1016/j.neuroimage.2019.116209>>\n"
    "       Heinzle et al. (2016) <<http://dx.doi.org/10.1016/j.neuroimage.2015.10.025>>\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char *f_input = NULL, *f_layer = NULL, *f_column = NULL, *f_ALF = NULL;
    char *f_out = NULL;
    string tag = "deveinDeconv";
    int ac;
    bool mode_linear = false;  // Default is linear as it requires least inputs
    bool mode_CBV = false;

    // Peak to tail ratio from Markuerkiaga et al. 2016 Fig. 5B at 7T.
    float lambda = 0.25;

    if (argc < 3) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layer_file\n");
                return 1;
            }
            f_layer = argv[ac];
        } else if (!strcmp(argv[ac], "-column_file")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -column_file\n");
                return 1;
            }
            mode_linear = false;
            f_column = argv[ac];
        } else if (!strcmp(argv[ac], "-ALF")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -ALF\n");
                return 1;
            }
            mode_linear = false;
            f_ALF = argv[ac];
        } else if (!strcmp(argv[ac], "-lambda")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -lambda\n");
                return 1;
            }
            lambda = atof(argv[ac]);  // No string copy, pointer assignment
        } else if (!strcmp(argv[ac], "-CBV")) {
            mode_CBV = true;
        } else if (!strcmp(argv[ac], "-linear")) {
            mode_linear = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            f_out = argv[ac];
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            f_input = argv[ac];
            f_out = argv[ac];  // Changes later if output is given.
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            f_out = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!f_input) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!f_layer) {
        fprintf(stderr, "** missing option '-layer_file'\n");
        return 1;
    }

    // Read inputs including data
    nifti_image* nii = nifti_image_read(f_input, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_input);
        return 2;
    }
    nifti_image* nii_layeri = nifti_image_read(f_layer, 1);
    if (!nii_layeri) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_layer);
        return 2;
    }

    log_welcome("LN2_DEVEIN");
    log_nifti_descriptives(nii);
    log_nifti_descriptives(nii_layeri);

    if (mode_linear) {
        tag = "deveinLinear";
        cout << "  Linear mode is selected!" << endl;
        cout << "    There will be no deconvolution, just layer-dependent depth scaling.\n" << endl;
    }
    if (mode_CBV) {
        tag = "deveinCBV";
        cout << "  CBV mode is selected!" << endl;
        cout << "    There will be no deconvolution, just layer-dependent CBV scaling.\n" << endl;
    }

    if (mode_linear && mode_CBV) {
        cout << "  You selected both: Linear and CBV scaling. So I am confused. " << endl;
        cout << "    Using only use linear scaling. \n" << endl;
         tag = "deveinLinear";
         mode_CBV = true;
    }

    // Get dimensions of input
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int size_z = nii->nz;
    const int size_time = nii->nt;
    const int nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii);
    float *nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* nii_layer = copy_nifti_as_float32(nii_layeri);
    float *nii_layer_data = static_cast<float*>(nii_layer->data);

    // Allocate new niftis
    nifti_image *nii_output = copy_nifti_as_float32(nii_input);
    float *nii_output_data = static_cast<float*>(nii_output->data);
    for (int i = 0; i < nr_voxels * size_time; ++i) {
        *(nii_output_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // Find input value ranges
    // ------------------------------------------------------------------------
    float layer_max = std::numeric_limits<float>::min();
    float layer_min = std::numeric_limits<float>::max();
    float input_max = std::numeric_limits<float>::min();
    float input_min = std::numeric_limits<float>::max();
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_layer_data + i) != 0) {
            if (*(nii_layer_data + i) > layer_max) {
                layer_max = *(nii_layer_data + i);
            }
            if (*(nii_layer_data + i) < layer_min) {
                layer_min = *(nii_layer_data + i);
            }
            if (*(nii_input_data + i) > input_max) {
                input_max = *(nii_input_data + i);
            }
            if (*(nii_input_data + i) < input_min) {
                input_min = *(nii_input_data + i);
            }
        }
    }

    // ------------------------------------------------------------------------
    // Find number of layers
    // ------------------------------------------------------------------------
    int nr_layers = 0;
    if (layer_max <= 1.0) {
        cout << "  !!! Normalized cortical depth metric is used." << endl;
        cout << "    Min. normalized cortical depth = " << layer_min << endl;
        cout << "    Max. normalized cortical depth = " << layer_max << endl;
    } else {
        for (int i = 0; i < nr_voxels; ++i) {
            if (*(nii_layer_data + i) > nr_layers) {
                nr_layers = *(nii_layer_data + i);
            }
        }
        cout << "  Number of discete layers = " << nr_layers << endl;
        cout << "    Min. layer value = " << layer_min << endl;
        cout << "    Max. layer value = " << layer_max << endl;
    }
    cout << endl;

    cout << "  Before devein:" << endl;
    cout << "    Min. value = " << input_min << endl;
    cout << "    Max. value = " << input_max << endl;
    cout << endl;

    // ========================================================================
    // Main chunk
    // ========================================================================
    if (mode_linear) {  // Handle linear case straightforwardly
        for (int t=0; t<size_time; ++t) {
            for (int i=0; i<nr_voxels; ++i) {
                int j = t * nr_voxels + i;
                if (*(nii_layer_data + i) > 0) {
                    float l;
                    if (layer_max <= 1.0) {  // Indicates metric file is used
                        l = 1. - *(nii_layer_data + i);
                    } else {
                        l =  1. - (*(nii_layer_data + i)-0.5)/layer_max ; // bug fix by ppxdm4 https://github.com/layerfMRI/LAYNII/issues/31
                    }

                    *(nii_output_data + j) = *(nii_input_data + j) * l;
                }
            }
        }
    } else {

        // Check additional inputs
        if (!f_column) {
            fprintf(stderr, "** missing option '-column_file'\n");
            return 1;
        }
        if (!f_ALF) {
            fprintf(stderr, "** missing option '-ALF'\n");
            return 1;
        }

        // Read additional inputs
        nifti_image* nii_columni = nifti_image_read(f_column, 1);
        if (!nii_columni) {
            fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_column);
            return 2;
        }
        nifti_image* nii_ALFi = nifti_image_read(f_ALF, 1);
        if (!nii_ALFi) {
            fprintf(stderr, "** failed to read NIfTI from '%s'\n", f_ALF);
            return 2;
        }

        log_nifti_descriptives(nii_columni);
        log_nifti_descriptives(nii_ALFi);

        // Prepare additional inputs
        nifti_image* nii_column = copy_nifti_as_int32(nii_columni);
        int32_t *nii_column_data = static_cast<int32_t*>(nii_column->data);
        nifti_image* nii_ALF = copy_nifti_as_float32(nii_ALFi);
        float* nii_ALF_data = static_cast<float*>(nii_ALF->data);

        // --------------------------------------------------------------------
        // Find number of columns
        // --------------------------------------------------------------------
        int nr_columns = 0;
        for (int i = 0; i < nr_voxels; ++i) {
            if (*(nii_column_data + i) > nr_columns) {
                nr_columns = *(nii_column_data + i);
            }
        }
        cout << "  Number of columns = " << nr_columns << endl;

        // --------------------------------------------------------------------
        // Making sure that every column voxel has a layer
        // --------------------------------------------------------------------
        for (int ivox = 0; ivox < nr_voxels; ++ivox) {
            if (*(nii_column_data + ivox) > 0 &&  *(nii_layer_data + ivox) == 0 ){
                cout << "  Some column voxels don't have layers in this file.\n";
                cout << "  It is likely that you provided wrong data.\n";
                cout << "  Though, the program will proceed." << endl;
                *(nii_column_data + ivox) = 0;
            }
        }

        // ====================================================================
        // Do deconvolution column by column. First, I allocate all.
        // ====================================================================
        // float vec1[nr_layers][size_time]
        // float vec2[nr_layers][size_time]
        // float vecALF[nr_layers];
        // int vec_nr_voxels[nr_layers];
        std::vector<std::vector<float>> vec1( nr_layers, std::vector<float>( size_time ) );
        std::vector<std::vector<float>> vec2( nr_layers, std::vector<float>( size_time ) );
        std::vector<float> vecALF(nr_layers);
        std::vector<int> vec_nr_voxels(nr_layers);

        for (int i = 0; i < nr_layers; ++i) {
            for (int t = 0; t < size_time; t++) {
                vec1[i][t] = 0.;
                vec2[i][t] = 0.;
            }
            vecALF[i] = 0.;
            vec_nr_voxels[i] = 0;
        }

        // ====================================================================
        // Big loop across columns
        // ====================================================================
        for (int icol = 1; icol <= nr_columns; ++icol) {

            // Reset vector
            for (int i = 0; i < nr_layers; ++i) {
                for (int t = 0; t < size_time; t++){
                    vec1[i][t] = 0.;
                    vec2[i][t] = 0.;
                }
                vecALF[i] = 0.;
                vec_nr_voxels[i] = 0;
            }

            // Fill vector of column #icol
            for (int ivox = 0; ivox < nr_voxels; ++ivox) {
                if (icol == *(nii_column_data + ivox)) {
                    int i = *(nii_layer_data + ivox) - 1;  // current layer

                    if (!mode_linear) {
                        vecALF[i] += *(nii_ALF_data + ivox);
                    }

                    vec_nr_voxels[i] += 1;

                    for (int t = 0; t < size_time; t++) {
                        vec1[i][t] += *(nii_input_data + t * nr_voxels + ivox);
                    }
                }
            }

            // Get mean of values within column vector
            for (int i = 0; i < nr_layers; ++i) {
                for (int t = 0; t < size_time; t++) {
                    vec1[i][t] /= (float)vec_nr_voxels[i];
                }
                vecALF[i] /= (float)vec_nr_voxels[i];
            }

            // ----------------------------------------------------------------
            // Do voxel deconvolution
            // ----------------------------------------------------------------

            // Normalize amplitude of low frequencies (ALF)
            float ALF_sum = 0 ;
            for (int i = 0; i < nr_layers; ++i) {
                if (vec_nr_voxels[i] > 0) {
                    ALF_sum += vecALF[i];
                }
            }
            for (int i = 0; i < nr_layers; ++i) {
                if (vec_nr_voxels[i] > 0) {
                    vecALF[i] /= ALF_sum;
                }
            }

            for (int t = 0; t < size_time; t++){

                for (int i = 0; i < nr_layers; ++i) {
                    if (mode_CBV) {  // Just CBV normalization
                        if (vec_nr_voxels[i] > 0) {
                            vec2[i][t] = vec1[i][t] / vecALF[i] * (float)nr_layers;
                        }
                    } else {  // Deconvolution
                        // Macrovascular contribution value, that needs to be
                        // subtracted from the current voxel.
                        float sum = 0;
                        for (int j = i-1; j >= 0; --j) {
                            // This is the deconvolution, It is weighted with CBV.
                            // Lambda is the inverse of peak to tail ratio
                            // from from Markuerkiaga et al. 2016 Fig. 5B at 7T.
                            if (vec_nr_voxels[j] > 0) {
                                sum += vec1[j][t] / (float)nr_layers / vecALF[j] * lambda;
                            }
                        }
                        if (vec_nr_voxels[i] > 0) {
                            vec2[i][t] = (vec1[i][t] - sum);
                        }
                    }
                }
            }

            // Fill file with the deconvolved values
            for (int ivox = 0; ivox < nr_voxels; ++ivox) {
                if (icol == *(nii_column_data + ivox)) {
                    int i = *(nii_layer_data + ivox) - 1;
                    for (int t = 0; t < size_time; t++) {
                        *(nii_output_data + t * nr_voxels + ivox ) = vec2[i][t];
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Find input value ranges
    // ------------------------------------------------------------------------
    input_max = std::numeric_limits<float>::min();
    input_min = std::numeric_limits<float>::max();
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_layer_data + i) != 0) {
            if (*(nii_output_data + i) > input_max) {
                input_max = *(nii_output_data + i);
            }
            if (*(nii_output_data + i) < input_min) {
                input_min = *(nii_output_data + i);
            }
        }
    }
    cout << "  After devein:" << endl;
    cout << "    Min. value = " << input_min << endl;
    cout << "    Max. value = " << input_max << endl;
    cout << endl;

    save_output_nifti(f_out, tag, nii_output, true);

    cout << "  Finished." << endl;
    return 0;
}
