#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_NOISE_KERNEL: Estimates functional point spread function by looking\n"
    "                 at shared sources of variance across different directions.\n"
    "\n"
    "Usage:\n"
    "    LN_NOISE_KERNEL -input timeseries.nii -kernel_size 11 \n"
    "\n"
    "test usage in the test_data folder: \n"
    "    ../LN_NOISE_KERNEL -input lo_Nulled_intemp.nii -kernel_size 7 \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -input       : Nifti (.nii) time series.\n"
    "    -kernel_size : (Optional) Use an odd positive integer (default 11).\n"
    "    -output      : (Optional) Output filename. Overwrites existing files.\n"
    "                   If not given, the prefix 'fPSF' is added \n"
    "\n"
    "Notes:\n"
    "    Some applications of this program are mentioned in this blog post: \n"
    "    <http://layerfmri.com/QA> \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ;
    char *fin = NULL;
    int ac;
    int kernel_size = 11; // This is the maximal number of layers. I don't know how to allocate it dynamically. this should be an odd number. That is smaller than half of the shortest matrix size to make sense
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
            kernel_size = atoi(argv[ac]);  // No string copy, pointer assignment
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


if (kernel_size%2==0) {
    cout << " you chose a kernel size of " << kernel_size << " even though I tols you to use an odd value... SHAME ON YOU" ;
    kernel_size = kernel_size -1 ;
    cout << "    I am using " << kernel_size << " instead" << endl;
}

    int kernel_vol = kernel_size * kernel_size* kernel_size ;
    double Nkernel[kernel_size][kernel_size][kernel_size] ;
    int Number_of_averages = 0;
    double Number_AVERAG[kernel_size][kernel_size][kernel_size] ;

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


for(int i = 0; i < kernel_size; i++) {
    for(int j = 0; j < kernel_size ; j++) {
        for(int k = 0; k < kernel_size ; k++) {
            Nkernel[i][j][k] = 0;
            Number_AVERAG[i][j][k] = 0;
        }
    }
}


    // Allocate new nifti
    nifti_image* nii_kernel = nifti_copy_nim_info(nii);
    nii_kernel->nt = 1;
    nii_kernel->nx = kernel_size;
    nii_kernel->ny = kernel_size;
    nii_kernel->nz = kernel_size;
    nii_kernel->nvox =  kernel_vol;
    nii_kernel->datatype = NIFTI_TYPE_FLOAT32;
    nii_kernel->nbyper = sizeof(float);
    nii_kernel->data = calloc(nii_kernel->nvox, nii_kernel->nbyper);
    float* nii_kernel_data = static_cast<float*>(nii_kernel->data);
    nii_kernel->scl_slope = 1; // to make sure that the units are given in Pearson correlations (-1...1)
    for(int ivoxel=0; ivoxel<kernel_vol; ++ivoxel) *(nii_kernel_data + ivoxel ) = 0.09;
    int knx = nii_kernel->nx;
    int knxy = nii_kernel->nx * nii_kernel->ny;
    int knxyz = nii_kernel->nx * nii_kernel->ny * nii_kernel->nz;


    // ========================================================================

cout << " Kernel size = " << kernel_size << endl;
cout << " Kernel size/2 = " << kernel_size/2 << endl;

//cout << " float = " << sizeof(float32_t) << endl;
//cout << " double = " << sizeof(double) << endl;

if ( size_x < kernel_size*2 || size_y < kernel_size*2 || size_z < kernel_size*2) {
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
        for(int kernely= -1*kernel_size/2; kernely<=kernel_size/2; ++kernely){
             for(int kernelx= -1*kernel_size/2; kernelx<=kernel_size/2; ++kernelx){
                 for(int kernelz= -1*kernel_size/2; kernelz<=kernel_size/2; ++kernelz){

                    vinc_x = ix + kernelx ;
                    vinc_y = iy + kernely ;
                    vinc_z = iz + kernelz ;
                    kern_ix = kernelx+kernel_size/2;
                    kern_iy = kernely+kernel_size/2;
                    kern_iz = kernelz+kernel_size/2;

                    if (vinc_x >= 0 && vinc_x < size_x && vinc_y >= 0 && vinc_y < size_y && vinc_z >= 0 && vinc_z < size_z) {

                        for(int it = 0 ; it < size_time  ; it++) {
                           vec2[it] = (double) *(nii_data  + nxyz *it +  nxy*vinc_z + nx*vinc_y + vinc_x) ;
                        }
                        dummy  = ren_correl(vec1, vec2, size_time) ;

                        //dummy = kernely ;
                        if (isfinite(dummy) && dummy != 0 ) {
                          //*(nii_kernel_data + knxy*kern_iz + knx*kern_iy + kern_ix) += dummy;
                           Nkernel[kern_iz][kern_iy][kern_ix] = Nkernel[kern_iz][kern_iy][kern_ix] +  dummy ;
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


        for(int kernely= 0; kernely<kernel_size; ++kernely){
             for(int kernelx= 0; kernelx<kernel_size; ++kernelx){
                 for(int kernelz= 0; kernelz<kernel_size; ++kernelz){

                    if (Number_AVERAG[kernelz][kernely][kernelx] > 0 ){
                       *(nii_kernel_data + knxy*kernelz + knx*kernely + kernelx) = (float) (Nkernel[kernelz][kernely][kernelx] / Number_AVERAG[kernelz][kernely][kernelx]); //dummy;
                    }
                    else *(nii_kernel_data + knxy*kernelz + knx*kernely + kernelx) =  0;
                 }
              }
          }


    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "fPSF", nii_kernel, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
