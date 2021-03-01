#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_MASK: This program is intended to help with voxels section for layer-profile extraction\n"
    "                    It generatees binary masks that span across all cortical depths.\n"
    "                    This is done to minimize circularity of only extracting layer signals.\n"
    "                    from voxels that exceed a significant detection threshold. \n"
    "                    Instead, if a voxel in a column is activated, \n"
    "                    the  signal of the entire colum can be extracted \n"
    "\n"
    "Usage:\n"
    "    LN2_MASK -scores activation.nii -columns columns.nii -max_thresh 2.3 -output column_mask.nii \n"
    "    ../LN2_MASK -scores activation.nii -columns columns.nii -mean_thresh 1.5 -output column_mask.nii \n"
    "    ../LN2_MASK -scores lo_BOLD_act.nii -columns lo_columns.nii -mean_thresh 1 -output mask.nii -abs \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -scores      : 3D Nifti image fucntional activation scores that will \n"
    "                   be transformed into layer dimensions.\n"
    "                   Those can be z-score maps of beta maps.\n"
    "                   The default is to exclusively consider at positive scores.\n"
    "    -columns     : A 3D nifti file that contains columns as intager masks.\n"
    "                   e.g. the output of LN2_COLUMNS or LN2_MULTILATERATE.\n"
    "    -max_thresh  : (Optional) Threshold of activation score. If any voxel in the column \n"
    "                   exceeds this value, the entire column in selected .\n"
    "                   This selection is use by default with value 1.0 \n"
    "    -mean_thresh : (Optional) Threshold of activation score. If the mean of the column  \n"
    "                   exceeds this value, the entire column in selected.\n"
    "                   if this parameter is used, the max_thresh option is ignored \n"
    "    -output      : (Optional) Output filename, including .nii or\n"
    "                   .nii.gz, and path if needed. Overwrites existing files.\n"
    "    -abs         : (Optional) if you want to also consider negative score values\n"
    "                   use this option.\n"
    "\n"
    "Notes:\n"
    "    - This refers to layerfMRI artifact here: \n"
    "    - Note that columns and scores are expected to be positive numbers \n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL;
    int ac;
    bool mode_max = true, mode_mean = false, mode_abs = false, use_outpath = false;
    float thresh = 1.0 ; 

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-scores")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -scores\n");
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
        } else if( !strcmp(argv[ac], "-mean_thresh") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -mean_thresh\n");
                return 1;
            }
            thresh = atof(argv[ac]);
            mode_mean = true; 
            mode_max = false; 
        } else if( !strcmp(argv[ac], "-max_thresh") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -max_thresh\n");
                return 1;
            }
            thresh = atof(argv[ac]);
            mode_max = true; 
        } else if( !strcmp(argv[ac], "-abs") ) {
            mode_abs = true; 
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            use_outpath = true;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-scores'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-columns'\n");
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

  
    if (mode_max == mode_mean) {
        cout << "There is something wrong, you need to select either mean thresholding or max thresholding. Not both." << endl;
        return 1; 
    }

    log_welcome("LN2_MASK");
    log_nifti_descriptives(nii1); //values
    log_nifti_descriptives(nii2); //columns
    
    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;

    const int nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];
    
    int nx = nii1->nx;
    int nxy = nii1->nx * nii1->ny;
    int nxyz = nii1->nx * nii1->ny * nii1->nz;


    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* columns = copy_nifti_as_int16(nii2);
    int16_t* columns_data = static_cast<int16_t*>(columns->data);
 
    // ========================================================================
    // making sure that there is nothing weird with the slope ]of the nii header 
    // ========================================================================

    if (nii_input->scl_slope == 0) {
        cout << "   It seems like the slope of the value file is ZERO" << endl;
        cout << "   I am setting it to 1 instead " << endl;
        nii_input->scl_slope = 1; 
    }



    // ========================================================================
    // Looking how many layers and columns we have and allocating the arrays accordingly
    // ========================================================================
    int nr_columns = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(columns_data + i) >= nr_columns){
            nr_columns = *(columns_data + i);
        }
    }
    cout << "    There are " << nr_columns<< " columns, total. " << endl << endl;


    double numb_voxels[nr_columns] ; 
    bool  thresh_exeed[nr_columns] ; 
    double mean_val[nr_columns] ; 
        for (int j = 0; j < nr_columns; j++) {
            mean_val   [j] = 0.; 
            numb_voxels[j] = 0.; 
            thresh_exeed[j] = false; 
        }

        // ========================================================================
        // considering necative activation too 
        // ========================================================================
    if (mode_abs) {
        for (int i = 0; i != nr_voxels; ++i) *(nii_input_data + i) = abs (*(nii_input_data + i));
    }
        // ========================================================================
        // Thresholding each column based on its maximally activated voxel
        // ========================================================================
    if (mode_max) {
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(nii_input_data + i) * nii_input->scl_slope  >= thresh){
                thresh_exeed [*(columns_data + i)-1] = true ;
            }
        }
    }

        
        // ========================================================================
        // Thresholding each column based on its mean activated signal within the colums
        // ========================================================================
    if (mode_mean) {
        
        // ========================================================================
        // average within columns and see if average is above threshold
        // ========================================================================
        for (int i = 0; i != nr_voxels; ++i) {
            if ( *(columns_data + i) > 0 ){
                mean_val   [ *(columns_data + i) -1 ] += *(nii_input_data + i) ;
                numb_voxels[ *(columns_data + i) -1 ] += 1; 
            }
        }
    
        for (int j = 0; j < nr_columns; j++) {
            if (numb_voxels[j] != 0){
                mean_val[j] /= (float)numb_voxels[j] ; 
                // cout << "layer " << i+1 << " and column " << j+1 << " has value " << mean_val[i][j] << endl;
            }
            
            if (mean_val[j] >= thresh){
                thresh_exeed [j] = true ; 
            }
        }
        
        
        
    }

    // ========================================================================
    // Prepare outputs
    nifti_image* mask = copy_nifti_as_int16(columns);
    int16_t* mask_data = static_cast<int16_t*>(mask->data);
    for (int voxi = 0; voxi <  nr_voxels; voxi++) *(mask_data + voxi) = 0.0 ; 



    // ========================================================================
    // setting all voxels in columns that have been selected
    // ========================================================================
    for (int i = 0; i < nr_voxels; ++i) {
        if ( (thresh_exeed [*(columns_data + i)-1] == true ) && *(columns_data + i) > 0 ){
           *(mask_data + i) =  1; 
        }
    }
    
 


    // ========================================================================
    // writing output
    // ========================================================================
    if (!use_outpath) fout = fin1;
    save_output_nifti(fout, "mask", mask, true, use_outpath);


  
    cout << "\n  Finished." << endl;
    return 0;
}
