#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_NOISE_KERNEL: \n"
    "         This program aims to estimate the functional point spread function\n"
    "         by looking at shared sources of variance across different directions.\n"
    "\n"
    "Usage:\n"
    "    LN_NOISE_KERNEL -input timeseries.nii -kernel_size 11 \n"
    "\n"
    "test usage in the test_data folder: \n"
    "    ../LN_NOISE_KERNEL -input lo_Nulled_intemp.nii -kernel_size 7 \n"
    "\n"
    "A few potential applications of this program are mentioned in this blog post: \n"
    "    http://layerfmri.com/QA \n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti (.nii) time series.\n"
    "    -kernal_size   : optional parameter for kernel size, use an odd positive intager (default 11)\n"
    "    -output        : (Optional) Custom output name. \n"
    "                     including the path, if you want to write it as specific locations \n"
    "                     including the file extension: nii or nii.gz \n"
    "                     This will overwrite excisting files with the same name \n"
    "                     if not given, the prefix fPSF is added \n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ; 
    char *fin = NULL;
    int ac;
    int kernal_size = 11; // This is the maximal number of layers. I don't know how to allocate it dynamically. this should be an odd number. That is smaller than half of the shortest matrix size to make sense
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-kernel_size")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -kernel_size\n");
                return 1;
            }
            kernal_size = atoi(argv[ac]);  // No string copy, pointer assignment
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        }
    }
    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }



    // Read input dataset
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_NOISE_KERNEL");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);


if (kernal_size%2==0) {
    cout << " you chose a kernel size of " << kernal_size << " even though I tols you to use an odd value... SHAME ON YOU" ; 
    kernal_size = kernal_size -1 ; 
    cout << "    I am using " << kernal_size << " instead" << endl;
}
    
    int kernal_vol = kernal_size * kernal_size* kernal_size ; 
    double Nkernal[kernal_size][kernal_size][kernal_size] ; 
    int Number_of_averages = 0; 
    double Number_AVERAG[kernal_size][kernal_size][kernal_size] ; 

/////////////////////////////////
////////allokate and se zero ////
/////////////////////////////////

double vec1[size_time] ;
double vec2[size_time] ;
double dummy = 0; 

    for(int timestep = 0; timestep < size_time ; timestep++) {
        vec1[timestep] = 0; 
        vec2[timestep] = 0; 
    }


for(int i = 0; i < kernal_size; i++) {
    for(int j = 0; j < kernal_size ; j++) {
        for(int k = 0; k < kernal_size ; k++) {
            Nkernal[i][j][k] = 0; 
            Number_AVERAG[i][j][k] = 0; 
        }
    }
}


    // Allocate new nifti
    nifti_image* nii_kernel = nifti_copy_nim_info(nii);
    nii_kernel->nt = 1;
    nii_kernel->nx = kernal_size;
    nii_kernel->ny = kernal_size;
    nii_kernel->nz = kernal_size;
    nii_kernel->nvox =  kernal_vol;
    nii_kernel->datatype = NIFTI_TYPE_FLOAT32;
    nii_kernel->nbyper = sizeof(float);
    nii_kernel->data = calloc(nii_kernel->nvox, nii_kernel->nbyper);
    float* nii_kernel_data = static_cast<float*>(nii_kernel->data);
    nii_kernel->scl_slope = 1; // to make sure that the units are given in Pearson correlations (-1...1) 
    for(int ivoxel=0; ivoxel<kernal_vol; ++ivoxel) *(nii_kernel_data + ivoxel ) = 0.09; 
    int knx = nii_kernel->nx;
    int knxy = nii_kernel->nx * nii_kernel->ny;
    int knxyz = nii_kernel->nx * nii_kernel->ny * nii_kernel->nz;


    // ========================================================================

cout << " Kernel size = " << kernal_size << endl; 
cout << " Kernel size/2 = " << kernal_size/2 << endl; 

//cout << " float = " << sizeof(float32_t) << endl; 
//cout << " double = " << sizeof(double) << endl; 

if ( size_x < kernal_size*2 || size_y < kernal_size*2 || size_z < kernal_size*2) {
    cout << "####################################################" << endl; 
    cout << "#### WARNING your Kernel might be too big ##########" << endl; 
    cout << "####################################################" << endl; 
} 

int vinc_x = 0 , vinc_y = 0 , vinc_z = 0; 
int kern_ix, kern_iy, kern_iz ; 


// four time estimate 
int all_loops = size_y * size_x ; 
int loop_counter = 0; 

for(int iy=0; iy<size_y; ++iy){
    for(int ix=0; ix<size_x; ++ix){
       cout << "\r"<<  loop_counter*100/all_loops  << "    % done "  << flush ; 
       loop_counter++ ;  
       for(int iz=0; iz<size_z; ++iz){           
               
          for(int it = 0 ; it < size_time  ; it++) {
               //cout  << "  it=" << it << "   ix="  <<ix << "   iy="  <<iy << "   iz="  <<iz<<   endl; 
               vec1[it] =  (double)*(nii_data  + nxyz *it +  nxy*iz + nx*iy + ix) ;


          }

          // going trhough vincinity of every voxel 
        for(int kernaly= -1*kernal_size/2; kernaly<=kernal_size/2; ++kernaly){
             for(int kernalx= -1*kernal_size/2; kernalx<=kernal_size/2; ++kernalx){
                 for(int kernalz= -1*kernal_size/2; kernalz<=kernal_size/2; ++kernalz){
                     
                    vinc_x = ix + kernalx ; 
                    vinc_y = iy + kernaly ; 
                    vinc_z = iz + kernalz ; 
                    kern_ix = kernalx+kernal_size/2; 
                    kern_iy = kernaly+kernal_size/2; 
                    kern_iz = kernalz+kernal_size/2; 

                    if (vinc_x >= 0 && vinc_x < size_x && vinc_y >= 0 && vinc_y < size_y && vinc_z >= 0 && vinc_z < size_z) {

                        for(int it = 0 ; it < size_time  ; it++) {
                           vec2[it] = (double) *(nii_data  + nxyz *it +  nxy*vinc_z + nx*vinc_y + vinc_x) ;
                        }
                        dummy  = ren_correl(vec1, vec2, size_time) ; 

                        //dummy = kernaly ; 
                        if (isfinite(dummy) && dummy != 0 ) {
                          //*(nii_kernel_data + knxy*kern_iz + knx*kern_iy + kern_ix) += dummy; 
                           Nkernal[kern_iz][kern_iy][kern_ix] = Nkernal[kern_iz][kern_iy][kern_ix] +  dummy ; 
                           Number_AVERAG[kern_iz][kern_iy][kern_ix]++ ; 

                        }
                    }
                     
                 }
              }
          }
         
         
         
       }
    }
}

  
cout << endl; 


        for(int kernaly= 0; kernaly<kernal_size; ++kernaly){
             for(int kernalx= 0; kernalx<kernal_size; ++kernalx){
                 for(int kernalz= 0; kernalz<kernal_size; ++kernalz){
                      
                    if (Number_AVERAG[kernalz][kernaly][kernalx] > 0 ){
                       *(nii_kernel_data + knxy*kernalz + knx*kernaly + kernalx) = (float) (Nkernal[kernalz][kernaly][kernalx] / Number_AVERAG[kernalz][kernaly][kernalx]); //dummy; 
                    }
                    else *(nii_kernel_data + knxy*kernalz + knx*kernaly + kernalx) =  0; 
                 }
              }
          }             


    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "fPSF", nii_kernel, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
