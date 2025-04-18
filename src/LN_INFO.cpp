#include "../dep/laynii_lib.h"

using namespace std;

void plotgray(double val, double mean, double stdev, bool inv );


int show_help(void) {
    printf(
    "LN_INFO: Displays the basic features of the nii image and plots some relevant header\n"
    "         information and attempts to get a terminal view of the image.\n"
    "\n"
    "Usage:\n"
    "    LN_INFO -input input_example.nii  \n"
    "    ../LN_INFO -input sc_UNI.nii -sub 5  \n" 
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
    cout << endl; 
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac; 
    int subs = 1;
    float std_val;
    bool NoPlotting = false;
    bool inv = false;  

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
            cout << "Not viewing the content of the BRIKS"  << endl;
        } else if (!strcmp(argv[ac], "-inv")) {
            inv = true;
            cout << "Using inverse colors"  << endl;
        } else if( ! strcmp(argv[ac], "-sub") ) {
          if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -sub\n");
            return 1;
           }
           subs = atof(argv[ac]);  // no string copy, just pointer assignment 
           cout << " Plotting the image by a factor of " << subs <<  " smaller" << endl;
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
    nifti_image* nii = nifti_image_read(fin, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_INFO");

    log_nifti_descriptives(nii);
    cout << "    Datatype = "  << nifti_datatype_string(nii->datatype) << endl; 
    cout << "    Slope = " << nii->scl_slope << endl; 
    cout << "    Intersept = " << nii->scl_inter << endl; 
    cout << "    Intent = code:" << nii->intent_code << " string:" << nifti_intent_string(nii->intent_code) << endl;  

    int64_t sizeSlice = nii->nz; 
    int64_t sizePhase = nii->ny; 
    int64_t sizeRead = nii->nx; 
    int64_t nrep = nii->nt; 
    int64_t nx = nii->nx;
    int64_t nxy = nii->nx * nii->ny;
    int64_t nxyz = nii->nx * nii->ny * nii->nz;

    // ========================================================================
    // Allocating new nifti
    nifti_image* nii_new = copy_nifti_as_float32(nii);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    // ========================================================================
    // Signal characteristics. 
    cout << endl << endl <<  "    BRIK value characeristics" << endl;  

    // Max value
    double max_val = -2147483648; 
    int64_t max_x = -1; 
    int64_t max_y = -1; 
    int64_t max_z = -1; 
    int64_t max_t = -1; 
    for(int64_t it=0; it < nrep; ++it){  
        for(int64_t iz=0; iz < sizeSlice; ++iz){  
            for(int64_t iy=0; iy < sizePhase; ++iy){
                for(int64_t ix=0; ix < sizeRead-0; ++ix){
                    if ( *(nii_new_data + nxyz*it + nxy*iz + nx*iy + ix) > max_val ) {
                        max_val = *(nii_new_data + nxyz *it + nxy*iz + nx*iy + ix);
                        max_x = ix; 
                        max_y = iy;
                        max_z = iz;
                        max_t = it;
                    }
                }
            }
        }
    }
    cout << "    Maximal value is " << max_val << " at (t=" << max_t << ",x=" << max_x << ",y=" << max_y << ",z=" << max_z << ")" << endl;

    // Min value
    double min_val = 2147483648; 
    int64_t min_x = -1; 
    int64_t min_y = -1; 
    int64_t min_z = -1; 
    int64_t min_t = -1; 
    for(int64_t it=0; it < nrep; ++it) {  
        for(int64_t iz=0; iz < sizeSlice; ++iz) {  
            for(int64_t iy=0; iy < sizePhase; ++iy) {
                for(int64_t ix=0; ix < sizeRead-0; ++ix) {
                    if ( *(nii_new_data + nxyz *it + nxy*iz + nx*iy + ix) < min_val ) {
                        min_val = *(nii_new_data + nxyz *it + nxy*iz + nx*iy+ ix);
                        min_x = ix; 
                        min_y = iy;
                        min_z = iz;
                        min_t = it;
                    }
                }
            }
        }
    }
    cout << "    Minimal value is "  << min_val << " at (t="  << min_t << ",x=" << min_x<< ",y=" << min_y<< ",z=" << min_z << ")" << endl;

    // Average value
    std::vector<double> vec1(nxyz/100);
    for (int64_t iv = 0; iv < nxyz/100; ++iv) {
        vec1[iv] =  static_cast<double>(*(nii_new_data + iv*100)); 
    }
    double mean_val = ren_average(vec1.data(), nxyz/100); 
    double stdev_val = ren_stdev(vec1.data(), nxyz/100); 
    cout << "    Average value is "  << mean_val << " and STEDV across space and time is " << stdev_val << endl;  
    cout << endl << endl << "    Attempt of plotting in terminal" << endl;  
 
    if (!NoPlotting) {    
        if (subs<1) subs = 1;
            cout << "Press enter to go to the next slice" << endl; 
            cout << "Press ctr+c to exit" << endl; 
            getchar();

        for(int64_t iz=0; iz < sizeSlice; iz = iz+subs) { 
            printf("\033[2J");
            printf("\033[%d;%dH", 0, 0);
            for(int64_t iy=0; iy < sizePhase; iy = iy + 2*subs ){
                for(int64_t ix=0; ix < sizeRead; ix = ix + subs ){
                plotgray( *(nii_new_data + nxyz * 0 + nxy*iz + nx*iy + ix), mean_val, stdev_val, inv);         
            }
            cout  << "\n";
          }
          cout << "Press enter to go to the next slice" << endl; 
          cout << "Press ctr+c to exit" << endl; 
          getchar();
        }
    }
    return 0;
}

void plotgray(double val, double mean, double stdev, bool inv){
    double val_n = (val-mean) / (stdev); 
    // if       (val_n > 0.6745                     ) cout << " ";  
    // // the value of 0.6745 refers to the z-score that contains 25% of the voxels for Gausssian distributions
    // else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "░"; 
    // else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "▒";
    // else if  (val_n                     < -0.6745) cout << "▓"; 
    // else (val_n >= 0.0  && val_n < 0.25 ) cout << "▓"; 

    if (inv){
        if       (val_n > 0.6745                     ) cout << "@";  
        else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "*"; 
        else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "-";
        else if  (val_n                     < -0.6745) cout << " "; 
    }

    if (!inv) {
        if       (val_n > 0.6745                     ) cout << " ";  
        else if  (val_n >= 0.      && val_n < 0.6745 ) cout << "░"; 
        else if  (val_n >= -0.6745 && val_n < 0.     ) cout << "▒";
        else if  (val_n                     < -0.6745) cout << "▓";
    }
}
