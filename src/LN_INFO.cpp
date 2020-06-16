
#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_INFO: displays the basic features of the nii image.\n"
    "         it plots some relevant header information \n"
    "         and attempts to git a terminal view of the image \n"
    "\n"
    "Usage:\n"
    "    LN_INFO -input input_example.nii  \n"
    "    ../LN_INFO -input lo_T1EPI.nii \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Specify input dataset.\n"

    "\n"
    "\n");
    cout << endl ; 
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac;
    float std_val;
    
    

    // Process user options
    if (argc < 2) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
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
    nifti_image* nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_INFO");
    cout << "    File name: " << nii_input->fname << endl;
    cout << "    Image details: " << nii_input->nz << " Z | " << nii_input->nx << " X | " << nii_input->ny << " Y | " << nii_input->nt << " T " << endl;
    cout << "    Voxel size = " << nii_input->pixdim[1] << " x " << nii_input->pixdim[2] << " x " << nii_input->pixdim[3] << endl;
    cout << "    Datatype: code =" << nii_input->datatype << ", string="  << nifti_datatype_string(nii_input->datatype) << endl; 
    cout << "    nii slope: "  << nii_input->scl_slope << endl; 
    cout << "    nii intersept : "  << nii_input->scl_inter << endl; 
    cout << "    nii intent: code="  << nii_input->intent_code << ", string="  << nifti_intent_string(nii_input->intent_code ) << endl;  

    
    // ========================================================================
    // Allocating new nifti
    nifti_image* nii_new = copy_nifti_as_float32(nii_input);
    float* nii_new_data = static_cast<float*>(nii_new->data);
    // ========================================================================
   int sizeSlice = nii_input->nz ; 
   int sizePhase = nii_input->ny ; 
   int sizeRead = nii_input->nx ; 
   int nrep =  nii_input->nt; 
   int nx =  nii_input->nx;
   int nxy = nii_input->nx * nii_input->ny;
   int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // signal characteristics. 


    cout << endl<< endl<<  "    BRIK value characeristics" << endl;  

    // max value
    double max_val = -2147483648 ; 
    int max_x = -1 ; 
    int max_y = -1 ; 
    int max_z = -1 ; 
    int max_t = -1 ; 
    
    for(int it=0; it<nrep; ++it){  
      for(int iz=0; iz<sizeSlice; ++iz){  
        for(int iy=0; iy<sizePhase; ++iy){
          for(int ix=0; ix<sizeRead-0; ++ix){
            if (*(nii_new_data    + nxyz *it + nxy*iz + nx*iy  + ix  ) > max_val ) {
                    max_val = *(nii_new_data    + nxyz *it + nxy*iz + nx*iy  + ix  ) ;
                    max_x = ix ; 
                    max_y = iy ;
                    max_z = iz ;
                    max_t = it ;
            }  
          }
        }
      }
    }
    cout << "    Maximal value is "  << max_val << " at (t="  << max_t<< ",x="<< max_x<< ",y="<< max_y<< ",z="<< max_z<< ")" << endl;  


    // min value
    double min_val = 2147483648 ; 
    int min_x = -1 ; 
    int min_y = -1 ; 
    int min_z = -1 ; 
    int min_t = -1 ; 
    
    for(int it=0; it<nrep; ++it){  
      for(int iz=0; iz<sizeSlice; ++iz){  
        for(int iy=0; iy<sizePhase; ++iy){
          for(int ix=0; ix<sizeRead-0; ++ix){
            if (*(nii_new_data    + nxyz *it + nxy*iz + nx*iy  + ix  ) < min_val ) {
                    min_val = *(nii_new_data    + nxyz *it + nxy*iz + nx*iy  + ix  ) ;
                    min_x = ix ; 
                    min_y = iy ;
                    min_z = iz ;
                    min_t = it ;
            }  
          }
        }
      }
    }
    cout << "    Minimal value is "  << min_val << " at (t="  << min_t<< ",x="<< min_x<< ",y="<< min_y<< ",z="<< min_z<< ")" << endl;  

    // average value
    double vec1[nxyz*nrep];
    for (int iv = 0; iv < nxyz*nrep; ++iv) {
        vec1[iv] =  static_cast<double>(*(nii_new_data + iv)) ; 
    }
    cout << "    Average value is "  << ren_average(vec1, nxyz*nrep) << " and STEDV across space and time is " << ren_stdev(vec1, nxyz*nrep)<< endl;  


    cout << endl<< endl<<  "    Attempt of plotting in terminal" << endl;  
 

    cout << "\n    Axial first TR" << endl;


      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
          *(nii_new_data    + nxyz *0 + nxy*sizeSlice/2 + nx*ix  + iy  )   ; 

        }
      }

   
    cout << "    Finished." << endl;
    return 0;
}

