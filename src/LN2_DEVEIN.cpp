
// TODO(Faruk): Seems there might be an issue with the gaussian kernel's
// symmetry and size corresponding to what is written in CLI

// To do make the vincinity direction specific vinc_x, vinc_y, vinc_z


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_DEVEIN : Program that tries to remove the macrovascular component in BOLD.\n"
    "\n"
    "    This program tries to estimate the micrcovascular component in layer-fMRI GE-BOLD \n"
    "    It does so by using an estimate of the local macrovascular blood volume ALF \n"
    "\n"
    "\n"
    "Usage:\n"
    "    LN2_DEVEIN -layer_file layers.nii -column_file columns.nii -input input time_series.nii -ALF ALF.nii  -FWHM 1\n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -layer_file    : Nifti (.nii) file that contains layers. \n"
    "    -column_file   : Nifti (.nii) file that should be smooth.  \n"
    "    -lambda        : Optional for the peak to tail ratio \n"
    "                     Default is 0.25 from Makuerikiaga Fig. 5B at 7T.\n"
    "                     If you assume that your CBF is exeptionally low, use 0.3, \n"
    "                     If you assume that your CBV is exeptionally high, use 0.2 \n"
    "                     If you use 3T instead of 7T, use 20 percent larger values \n"
    "                     Larger values will result in stronger deconvolution \n"
    "                     These variations will not affect the resulting profiles a lot\n"
    "                     So I would recommend to not touch this parameter to begin with\n"
    "    -linear        : Optional flag for linear layer scaling \n"
    "                     This means that there will be no deconvolution, just layer-dependent scaling\n"
    "    -CBV           : Optional flag for CBV scaling \n"
    "                     This means that there will be no deconvolution, just scaling of the baseline venous CBV\n"
    "    -input         : BOLD file that should be corrected from macro-vascular contaminations\n"
    "                     This can be a time series or an activity map (not z-scores though). \n"
    "    -ALF           : File with estiates of Amplitude of low frequencies as a correlate ot venous CBV  \n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* f_input = NULL, *f_layer = NULL, *f_column = NULL, *f_ALF = NULL;
    int ac, do_masking = 0, sulctouch = 0;
    float FWHM_val = 0;
    bool linear = false; 
    bool CBV = false; 
    float lambda = 0.25 ; // this is the peak to tail ratio from Makuerikiaga Fig. 5B at 7T
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
            f_column = argv[ac];
        } else if (!strcmp(argv[ac], "-ALF")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -ALF\n");
                return 1;
            }
            f_ALF = argv[ac];
        } else if (!strcmp(argv[ac], "-lambda")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -lambda\n");
                return 1;
            }
            lambda = atof(argv[ac]);  // No string copy, pointer assignment
        } else if (!strcmp(argv[ac], "-linear")) {
            linear = true ; 
            cout << " there will be no deconvolution, just layer-dependent depth scaling" << endl; 
        } else if (!strcmp(argv[ac], "-CBV")) {
            CBV = true ; 
            cout << " there will be no deconvolution, just layer-dependent CBV scaling" << endl; 
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            f_input = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

   if (linear && CBV) {
        fprintf(stderr, " I don't know what you want me to do linear or CBV? Both is not possible at the same time");
        return 1;
    }

    if (!f_input) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!f_layer) {
        fprintf(stderr, "** missing option '-layer_file'\n");
        return 1;
    }
    if (!f_column) {
        fprintf(stderr, "** missing option '-column_file'\n");
        return 1;
    }
    if (!f_ALF) {
        fprintf(stderr, "** missing option '-ALF'\n");
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

    log_welcome("LN2_DEVEIN");
    log_nifti_descriptives(nii);
    log_nifti_descriptives(nii_layeri);
    log_nifti_descriptives(nii_columni);
    log_nifti_descriptives(nii_ALFi);

    // Get dimensions of input
    const int size_z = nii->nz;
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int nx = nii->nx;
    const int nxy = nii->nx * nii->ny;
    const int nr_voxels = size_z * size_y * size_x;
    const int nrep = nii->nt;
    const int nr_voxt = nr_voxels * nrep;
    const float dX = nii->pixdim[1];
    const float dY = nii->pixdim[2];
    const float dZ = nii->pixdim[3];

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii);
    float *nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* nii_layer = copy_nifti_as_int32(nii_layeri);
    int32_t *nii_layer_data = static_cast<int32_t*>(nii_layer->data);
    nifti_image* nii_column = copy_nifti_as_int32(nii_columni);
    int32_t *nii_column_data = static_cast<int32_t*>(nii_column->data);
    nifti_image* nii_ALF = copy_nifti_as_float32(nii_ALFi);
    float *nii_ALF_data = static_cast<float*>(nii_ALF->data);

    // Allocate new niftis
    nifti_image *nii_decov = copy_nifti_as_float32(nii_input);
    float *nii_decov_data = static_cast<float*>(nii_decov->data);


    // Zero new images
    for (int i = 0; i < nr_voxt; ++i) *(nii_decov_data + i) = 0;


    ///////////////////////////
    // Find number of layers //
    ///////////////////////////
    int nr_layers = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_layer_data + i) > nr_layers)  nr_layers = *(nii_layer_data + i);
    }
    cout << "  There are " << nr_layers << " layers to deconvolve." << endl;
    
    ///////////////////////////
    // Find number of columns //
    ///////////////////////////
    int nr_columns = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_column_data + i) > nr_columns)  nr_columns = *(nii_column_data + i);
    }
    cout << "  There are " << nr_columns << " columns to do the deconvolution within." << endl;
    
    
    /////////////////////////////////////////////////////////////
    //   I will do the deconvolution column by column ///////////
    //   first, I allocate all, I need to do the deconvolution //
    /////////////////////////////////////////////////////////////

    float vec1[nr_layers][nrep], vec2[nr_layers][nrep], vecALF[nr_layers];
    int nx_voxls[nr_layers];
    for (int i = 0; i < nr_layers; ++i) {
        for (int timestep = 0; timestep < nrep; timestep++){
             vec1[i][timestep] = 0.;  
             vec2[i][timestep] = 0.; 
         }
        vecALF[i] = 0.;
        nx_voxls[i] = 0; 
    }
    int cur_layer = 0; 
    float cur_ALFmean = 0.; 
    float sum = 0. ; // value of macrovascular contribution, that needs to be subtracted from the current voxel

    //////////////////////////////////////////////////////
    // Making sure that every column voxel has a layer  //
    //////////////////////////////////////////////////////
    for (int ivox = 0; ivox < nr_voxels; ++ivox) {
        if (*(nii_column_data + ivox) > 0 &&  *(nii_layer_data + ivox) == 0 ){
            cout << " You gave me a file where some column voxels don't have layers "<< endl;
            cout << " This meand that you likely gave me the wrong data" << endl; 
            cout << " Though, I will try to work with is anyway" << endl; 
            *(nii_column_data + ivox) = 0;
        }
    }



    ////////////////////////////
    // Big loop across colums //
    ////////////////////////////
    for (int icol = 1; icol <= nr_columns; ++icol) {
        
        // resetting vector
        for (int i = 0; i < nr_layers; ++i) {
            for (int timestep = 0; timestep < nrep; timestep++){
                 vec1[i][timestep] = 0.;  
                 vec2[i][timestep] = 0.; 
             }
            vecALF[i] = 0.;
            nx_voxls[i] = 0; 
         }
        
       // filling vector of column #icol
       for (int ivox = 0; ivox < nr_voxels; ++ivox) {
           if (icol == *(nii_column_data + ivox)) {
               cur_layer = *(nii_layer_data + ivox)-1 ; 
               vecALF[cur_layer] = vecALF[cur_layer] +  *(nii_ALF_data + ivox) ;
               nx_voxls[cur_layer]++ ; 
               for (int timestep = 0; timestep < nrep; timestep++) vec1[cur_layer][timestep] = vec1[cur_layer][timestep] + *(nii_input_data + timestep*nr_voxels + ivox ) ;
           }
        } 
        
        // getting mean of falues within column vector 
        for (int i = 0; i < nr_layers; ++i) {
                for (int timestep = 0; timestep < nrep; timestep++) vec1[i][timestep] = vec1[i][timestep]/(float)nx_voxls[i] ;  
                vecALF[i] = vecALF[i]/(float)nx_voxls[i];
        }
        
        
        ///////////////////////////////
        // doing voxel deconvolution
        ////////////////////////////
        
        // before output for debugging
        cout <<  "column " << icol << "    of " << nr_columns << " Nvoxles = " ;
        for (int i = 0; i < nr_layers; ++i) cout << nx_voxls[i] << "  " ;
        cout << endl;  
        
        cout <<  "column " << icol << "    of " << nr_columns << " ALF = " ;
        for (int i = 0; i < nr_layers; ++i) cout << vecALF[i] << "  " ;
        cout << endl;  
        
        cout <<  "column " << icol << "    of " << nr_columns << " vec1 = " ;
        for (int i = 0; i < nr_layers; ++i) cout << vec1[i][0] << "  " ; 
        cout << endl;  

               //doing normalization of ALF
               cur_ALFmean = 0 ; 
        for (int i = 0; i < nr_layers; ++i) {
            if (nx_voxls[i] > 0 ) cur_ALFmean = cur_ALFmean + vecALF[i];   
         }
        for (int i = 0; i < nr_layers; ++i) {
            if (nx_voxls[i] > 0 ) vecALF[i] = vecALF[i] /cur_ALFmean ;   
         }
         
        cout <<  "column " << icol << "    of " << nr_columns << " ALF normaliced = " ;
        for (int i = 0; i < nr_layers; ++i) cout << vecALF[i] << "  " ;
        cout << endl;  
        
        for (int timestep = 0; timestep < nrep; timestep++){
                           //inicializing vec 2
            //for (int i = 0; i < nr_layers; ++i) vec2[i][timestep] = vec1[i][timestep];   
    
                          //getting whieghted sum
            for (int i = 0; i < nr_layers; ++i) {
                sum = 0;
                for (int j = i-1; j >= 0; --j) {
                    // this is the actual deconvolution, It is whieghted with CBV. 
                    // the lambda is the inverse of the peak to tail ratio from Markuerikiaga Fig. 5B (7Tesla)
                    // 
                    if (nx_voxls[j] > 0 ) sum = sum + vec1[j][timestep]/(float)nr_layers/vecALF[j]*lambda  ; 
                }
                if (nx_voxls[i] > 0 ) vec2[i][timestep] = (vec1[i][timestep] - sum) ;  
                
                // if just linear correction
                if (nx_voxls[i] > 0 && linear) vec2[i][timestep] = vec1[i][timestep]/(float)(i+1)*nr_layers ;
                
                // if just CBV normalication
                if (nx_voxls[i] > 0 && CBV ) vec2[i][timestep] = vec1[i][timestep]/vecALF[i]*(float)nr_layers ;
            }
         }
         
         // after output for debugging
        cout <<  "column " << icol << "    of " << nr_columns << " vec2 = " ;
        for (int i = 0; i < nr_layers; ++i) cout << vec2[i][0] << "  " ; 
        cout << endl;  
  
        
        // filling the file with the deconvolved values 
        for (int ivox = 0; ivox < nr_voxels; ++ivox) {
           if (icol == *(nii_column_data + ivox)) {
               cur_layer = *(nii_layer_data + ivox)-1 ; 
               for (int timestep = 0; timestep < nrep; timestep++) *(nii_decov_data + timestep*nr_voxels + ivox ) = vec2[cur_layer][timestep] ;
           }
        } 
        
    //cout << "\r" << "column " << icol << "    of " << nr_columns << "   " <<  flush ;  
    cout << endl << endl; 
    }// loop acrtoss columns closed
    cout << endl; 
    


    save_output_nifti(f_input, "deconvolved", nii_decov, true);

    cout << "  Finished." << endl;
    return 0;
}
