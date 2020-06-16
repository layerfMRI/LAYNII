
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
   int sizePhase = nii_input->nx ; 
   int sizeRead = nii_input->ny ; 
   int nrep =  nii_input->nt; 
   int nx =  nii_input->nx;
   int nxy = nii_input->nx * nii_input->ny;
   int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // obtain good value range for plotting. 
    
    int nr_voxels = nii_input->nvox;
    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_new_data + i) += 1 ; 
    }

    cout << "\n    axial first TR" << endl;


for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
          *(nii_new_data    + nxyz *0 + nxy*iz + nx*ix  + iy  ) =  1.  ; 

        }
      }
    }

   
    cout << "    Finished." << endl;
    return 0;
}

