
// to be compiled without make file as follows : 
//  c++ -std=c++11 -DHAVE_ZLIB -o LN_INFO src/LN_INFO.cpp dep/nifti2_io.cpp dep/znzlib.cpp  -I  -lm -lz
// note the two dependencies of nifti1.h, nifti2_io.cpp, nifti2_io.h, nifti2.h, README.md, znzlib.cpp , znzlib.h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "../dep/nifti2_io.h"

using namespace std;

nifti_image* copy_nifti_as_float32(nifti_image* nii);
void plotgray(double val, double mean, double stdev, bool inv );
double ren_average(double arr[], int size);
double ren_stdev(double arr[], int size);

int show_help(void) {
    printf(
    "LN_INFO: displays the basic features of the nii image.\n"
    "         it plots some relevant header information \n"
    "         and attempts to git a terminal view of the image \n"
    "\n"
    "Usage:\n"
    "    compile this program with the following command: \n"
    "        c++  -DHAVE_ZLIB -o LN_INFO LN_INFO.cpp nifti2_io.cpp znzlib.cpp  -I  -lm -lz \n"
    "\n"
    "    Execute the program with the following command: \n" 
    "    LN_INFO -input input_example.nii  \n"
    "    ../LN_INFO -input lo_T1EPI.nii -sub 3  \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Specify input dataset.\n"
    "    -NoPlot : (Optional) In case you do not want to plot the content of the BRIKS.\n"
    "    -sub    : (Optional) subsample plotting to make it smaller.\n"
    "              the number given after -sub is the factor of voxels to skip \n" 
    "    -inv    : (Optional) invert color scale for black terminal.\n"
    "\n"
    "\n");
    cout << endl ; 
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac ; 
    int subs = 1. ;
    float std_val;
    bool NoPlotting = false ;
    bool inv = false ;  
    

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
        } else if (!strcmp(argv[ac], "-NoPlot")) {
            NoPlotting = true;
            cout << "I am not viewing the content of the BRIKS"  << endl;
        } else if (!strcmp(argv[ac], "-inv")) {
            inv = true;
            cout << "I am usinf inverse colors"  << endl;
        } else if( ! strcmp(argv[ac], "-sub") ) {
          if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -sub\n");
            return 1;
           }
           subs = atof(argv[ac]);  // no string copy, just pointer assignment 
           cout << " I will plor the imagien by a factor of " << subs <<  " smaller" << endl;
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
    //cout << "nxyz+"<< nxyz*nrep  << endl; 
    //nxyz = 100 ; 
    double vec1[nxyz/100];
    for (int iv = 0; iv < nxyz/100; ++iv) {
        vec1[iv] =  static_cast<double>(*(nii_new_data + iv*100)) ; 
    }
    double mean_val = ren_average(vec1, nxyz/100) ; 
    double stdev_val = ren_stdev(vec1, nxyz/100) ; 
    cout << "    Average value is "  << mean_val << " and STEDV across space and time is " << stdev_val << endl;  


    cout << endl<< endl<<  "    Attempt of plotting in terminal" << endl;  
 


  if (!NoPlotting) {
    //plotting
    
    if (subs<1) subs = 1;
    
      cout << "press enter to go to the next slice" << endl; 
      cout << "press ctr+c to exit" << endl; 
      getchar();
    
    for(int iz=0; iz<sizeSlice; iz= iz+subs){ 
      
      printf("\033[2J");
      printf("\033[%d;%dH", 0, 0);
       for(int iy=0; iy<sizePhase; iy = iy + 2*subs ){
          for(int ix=0; ix<sizeRead; ix = ix + subs ){
           //cout << ix <<  "  " << iy << "  " ;  
          plotgray( *(nii_new_data    + nxyz * 0 + nxy*iz + nx*iy  + ix  ), mean_val,  stdev_val, inv)    ; 
        
        }
        cout  << "\n"  ;
      }
      cout << "press enter to go to the next slice" << endl; 
      cout << "press ctr+c to exit" << endl; 
      getchar();
    }

  }
    


    
    return 0;
}

void plotgray(double val, double mean, double stdev, bool inv){
    double val_n = (val-mean) / (stdev) ; 
//    if       (val_n > 0.6745                     ) cout << " " ;  
//    // the value of 0.6745 refers to the z-score that contains 25% of the voxels for Gausssian distributions
//    else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "░" ; 
//    else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "▒" ;
//    else if  (val_n                     < -0.6745) cout << "▓" ; 
    ////    else (val_n >= 0.0  && val_n < 0.25 ) cout << "▓" ; 

    if (inv){
        if       (val_n > 0.6745                     ) cout << "@" ;  
        else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "*" ; 
        else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "-" ;
        else if  (val_n                     < -0.6745) cout << " " ; 
    }

    if (!inv) {
        if       (val_n > 0.6745                     ) cout << " " ;  
        else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "░" ; 
        else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "▒" ;
        else if  (val_n                     < -0.6745) cout << "▓" ;
    }

}


double ren_average(double arr[], int size) {
    int i;
    double sum = 0;
    double avg;

    for (i = 0; i < size; ++i) {
        sum += arr[i];
    }
    avg = double(sum) / size;
    return avg;
}

double ren_stdev(double arr[], int size) {
    int i;
    double sum = 0;
    double mean;

    mean = ren_average(arr, size);

    for (i = 0; i < size; ++i) {
        sum += (arr[i] - mean) * (arr[i] - mean) / ((double)size - 1);
    }
    return sqrt( sum);
}

nifti_image* copy_nifti_as_float32(nifti_image* nii) {
    ///////////////////////////////////////////////////////////////////////////
    // NOTE(Renzo): Fixing potential problems with different input datatypes //
    // here, I am loading them in their native datatype and cast them        //
    // to the datatype I like best.                                          //
    ///////////////////////////////////////////////////////////////////////////

    // NOTE(for future reference): Rick's comments:
    // nifti_copy_nim_info(). It will return with data == NULL.
    // If you need the data allocated, memory use would not change once you do
    // so. There is also nifti_make_new_nim()
    // nifti_image* nifti_make_new_nim(const int64_t dims[],
    //                                 int datatype, int data_fill)

    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    float* nii_new_data = static_cast<float*>(nii_new->data);
    int nr_voxels = nii_new->nvox;

    // NOTE(Faruk): See nifti1.h for notes on data types
    // ------------------------------------------------------------------------
    if (nii->datatype == 2) {  // NIFTI_TYPE_UINT8
        uint8_t* nii_data = static_cast<uint8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 512) {  // NIFTI_TYPE_UINT16
        uint16_t* nii_data = static_cast<uint16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 768) {  // NIFTI_TYPE_UINT32
        uint32_t* nii_data = static_cast<uint32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 1280) {  // NIFTI_TYPE_UINT64
        uint64_t* nii_data = static_cast<uint64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 256) {  // NIFTI_TYPE_INT8
        int8_t* nii_data = static_cast<int8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 4) {  // NIFTI_TYPE_INT16
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 8) {  // NIFTI_TYPE_INT32
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 1024) {  // NIFTI_TYPE_INT64
        int64_t* nii_data = static_cast<int64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 16) {  // NIFTI_TYPE_FLOAT32
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else if (nii->datatype == 64) {  // NIFTI_TYPE_FLOAT64
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    } else {
        cout << "Warning! Unrecognized nifti data type!" << endl;
    }

    // Replace nans with zeros
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_new_data + i)!= *(nii_new_data + i)) {
            *(nii_new_data + i) = 0;
        }
    }

    return nii_new;
}
