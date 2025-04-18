// WORK in PROGRESS
// TODO(Renzo): Think about what to do with time series data.
// E.g. carpet plot: https://github.com/layerfMRI/repository/tree/master/Layer_me

#include <fstream>
#include <iomanip>

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_PROFILE: Generates layer profiles from 3D nii file based on layer masks.\n"
    "             It averages all the signal intensities of each layer and write it\n"
    "             out as a 2d-plot. The output is a text file (table).\n"
    "               - Column 1 is the layer number.\n"
    "               - Column 2 is the mean signal in this layer.\n"
    "               - Column 3 is the STDEV of the signal variance across all voxels in this layer.\n"
    "               - Column 4 is the number of voxels per layer.\n"
    "\n"
    "Usage:\n"
    "    LN2_PROFILE -input activitymap.nii -layers layers.nii -plot \n"
    "    LN2_PROFILE -input activitymap.nii -layers layers.nii -plot -output layer_profile.txt \n"
    "    ../LN2_PROFILE -input sc_VASO_act.nii -layers sc_layers.nii -plot -debug \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -layers : Specify input dataset of layers\n"
    "              It is assumed that it conists of intager numbers of layers\n"
    "              It is assumed that deeper layers have small values \n"
    "              It is assumes that superficial layers have large values.\n"
    "    -input  : Specify input dataset of to extract the signal from.\n"
    "              This is usually an activation map.\n"
    "              This 3D nii file must have the same dimension as the layer file.\n"
    "    -mask   : (Optional) Specify a local mask, an ROI to pool the signal from.\n"
    "              This is usefull, if the layer input is larder than the ROI.\n"
    "              For many concentional pipelines this might be the output of LN2_MASK.\n"
    "              This 3D nii file must have the same dimension as the layer file.\n"
    "    -plot   : (Optional)\n"
    "              this option tries to plot the profile as ASKII art in the terminal \n"
    "              This option can be useful if you do not have a graphical plotting profile ready\n"
    "              E.g. on a remote server without X11 forwarding.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output basename.\n"
    "              Default is adding '_padded' as suffix \n"
    "\n"
    "Notes:\n"
    "    - The averaging is done across all voxels layers, independent of their value.\n"
    "    - If you only want to use average across a subset of layers, consider restricting the layer mask.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    uint16_t ac;
    nifti_image *nii1 = NULL;
    nifti_image *niil = NULL;
    nifti_image *niim = NULL;
    char *fin = NULL, *finl = NULL, *finm = NULL;
    char const *fout = "profile.txt";
    bool  mode_debug = false,  mode_plot = false;
    bool  use_outpath = false;
    bool  use_mask = false;

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
            fout = fin;
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layers\n");
                return 1;
            }
            finl = argv[ac];
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            use_mask = true; 
            finm = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-plot")) {
            mode_plot = true;
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finl) {
        fprintf(stderr, "** missing option '-layers'\n");
        return 1;
    }
    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }
    niil = nifti_image_read(finl, 1);
    if (!niil) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", finl);
        return 2;
    }
    
    if (use_mask == true) {
        niim = nifti_image_read(finm, 1);
        if (!niim) {
            fprintf(stderr, "** failed to read NIfTI from '%s'\n", finm);
            return 2;
        }
    }
    
    log_welcome("LN2_PROFILE");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;
    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Load input
    // ========================================================================
    nifti_image* layers = copy_nifti_as_int16(niil);
    int16_t* layers_data = static_cast<int16_t*>(layers->data);
     
    nifti_image* act = copy_nifti_as_float32(nii1);
    float* act_data = static_cast<float*>(act->data);

    // ========================================================================
    // Make sure that there is nothing weird with the slope of the nii header
    // ========================================================================
    if(mode_debug) {
       cout << "   Layer file has slope " << niil->scl_slope  << endl;
       cout << "   Act  file has slope  " << nii1->scl_slope  << endl;
       if (use_mask == true)  cout << "   Mask file has slope  " << niim->scl_slope  << endl;
    }

    if (act->scl_slope == 0) {
        cout << "   There is something weird with the slope of the act file " << endl;
        cout << "   it seems to be ZERO, this doesn't make sense " << endl;
        cout << "   I am setting it to 1 instead " << endl;
        act->scl_slope = 1;
    }
    
    // ========================================================================
    // remove voxels outside mask, if ther is one specified.
    // ========================================================================
    if (use_mask == true ) {
	  nifti_image* mask = copy_nifti_as_int16(niim);
      int16_t* mask_data = static_cast<int16_t*>(mask->data);	
      
	  for (int j = 0; j != nr_voxels; ++j) {
            if (*(mask_data + j) == 0  ) {
                *(layers_data + j) = 0;
            }
        }
    }
    
    
    // ========================================================================
    // Look how many layers we have and allocating the arrays accordingly
    // ========================================================================
    int nr_layers = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(layers_data + i) >= nr_layers) {
            nr_layers = *(layers_data + i);
        }
    }
    cout << "    There are " << nr_layers<< " layers. " << endl << endl;

    // double numb_voxels[nr_layers];
    // double mean_layers[nr_layers];
    // double std_layers [nr_layers];
    std::vector<double> numb_voxels(nr_layers);
    std::vector<double> mean_layers(nr_layers);
    std::vector<double> std_layers (nr_layers);

    for (int i = 0; i < nr_layers; i++) {
        mean_layers[i] = 0.;
        std_layers[i] = 0.;
        numb_voxels[i] = 0.;
    }
	  
    // ========================================================================
    // Look how many voxels we have per layer
    // ========================================================================
    for(int i = 0; i < nr_layers; i++) {
        for (int j = 0; j != nr_voxels; ++j) {
            if (*(layers_data + j) == i+1 ) {
                numb_voxels[i] ++;
            }
        }
    }

    //-------------- finding layer with maximal number of voxels
    int max_layer_number = 0;
    int max_layer_number_layer = 0;
    for(int i = 0; i < nr_layers; i++) {
        if (numb_voxels[i] >= max_layer_number ) {
            max_layer_number =  numb_voxels[i];
            max_layer_number_layer = i;
        }
    }

    if(mode_debug) cout << "   Layer  " <<   max_layer_number_layer+1 << " has the most voxels: " <<  max_layer_number << endl;

    // ========================================================================
    // Go through layers MAIN loop
    // ========================================================================
    // double vec1[max_layer_number];
    std::vector<double> vec1(max_layer_number);
    int dummy_index = 0;

    for(int i = 0; i < nr_layers; i++) {
        for (int j = 0; j != nr_voxels; ++j) {
            if (*(layers_data + j) == i+1 ) {
                vec1[dummy_index] = *(act_data + j);
                dummy_index ++;
            }
        }
        dummy_index = 0;
        mean_layers[i] = ren_average(vec1.data(), numb_voxels[i])*act->scl_slope;
        std_layers[i]  = ren_stdev  (vec1.data(), numb_voxels[i])*act->scl_slope;
    }

    // ========================================================================
    // Write layer profiles to terminal
    // ========================================================================
    if (mode_debug) {
        for(int i = 0; i < nr_layers; i++) {
            cout << "In layer " << i+1 << " with a mean signal of "<<  mean_layers[i];
            cout <<  " +/-  " << std_layers[i] << " with  " <<  numb_voxels[i] <<  " are voxels " << endl;
        }
    }

    // ========================================================================
    // Write layer profiles to text file with the right file name
    // ========================================================================
    // Managing file name, path and extension
    string path_out;
    const string path = fout;

    if (use_outpath) {
        path_out = path;
    } else {
        // Parse path
        string dir, file, basename, ext, sep;
        auto pos1 = path.find_last_of('/');
        if (pos1 != string::npos) {  // For Unix
            sep = "/";
            dir = path.substr(0, pos1);
            file = path.substr(pos1 + 1);
        } else {  // For Windows
            pos1 = path.find_last_of('\\');
            if (pos1 != string::npos) {
                sep = "\\";
                dir = path.substr(0, pos1);
                file = path.substr(pos1 + 1);
            } else {  // Only the filename
                sep = "";
                dir = "";
                file = path;
            }
        }

        // Parse extension
        auto const pos2 = file.find_first_of('.');
        if (pos2 != string::npos) {
            basename = file.substr(0, pos2);
            ext = ".txt";
        } else {  // Determine default extension when no extension given
            basename = file;
            ext = ".txt";
        }

        // Prepare output path
        path_out = dir + sep + basename + "_" + "profile" + ext;
    }

    // Writing into file
    ofstream outf(path_out);
    if (!outf) {
        cout<<"error when opening the text file"<<endl;
    }

    cout<<"    writing to disk "  << path_out<<endl;
    for(int i = 0; i < nr_layers; i++) {
      outf << i+1 << "   "<<  mean_layers[i] <<  " " << std_layers[i] << "  " <<  numb_voxels[i] << endl;
     }
    outf.close();

    // ========================================================================
    // Plot in terminal, use ASCII to avoid issues with terminal types
    // ========================================================================
    if (mode_plot) {
        double max_val = -3.4028234664e+38;
        double min_val = 3.4028234664e+38;

        for (int i = 0; i < nr_layers; i++) {
            if (mean_layers[i] <= min_val) min_val = mean_layers[i];
            if (mean_layers[i] >= max_val) max_val = mean_layers[i];
        }

    cout.precision(2);
    cout << endl<< endl;

    // terminal width. of course this can be set automatically, but then it
    // will get dependencies of operating system
    int termwdth = 80;
    int termhght = 20;  // terminal height
    // int matrix[termwdth][termhght];
    std::vector<std::vector<double>> matrix( termwdth, std::vector<double>( termhght ) );

    for (int w =0; w < termwdth; w++) {
        for (int h =0; h < termhght;h++) {
            matrix [w][h] = 0;
        }
    }

    // filling matrix
    // padding for visually pleasing
    double max_valp = max_val + 1./(double)termwdth * (max_val-min_val);
    double min_valp = min_val - 1./(double)termwdth * (max_val-min_val);

    double x_lay = 0.;
    double y_val = 0.;

    for (int w =0; w < termwdth; w++) {
        for (int h =0; h < termhght;h++) {

            x_lay = (double) w / (double) termwdth * (double) nr_layers;
            y_val = (mean_layers[(int)x_lay] - min_valp) / (max_valp-min_valp) * termhght;

            if ( ( h - (int)y_val ) < 1 ) {
                matrix [w][h] = 1;
            }
        }
    }

    // top bar // two lines
    cout << "       +-";
    for (int w =0; w < termwdth; w++) cout << "-";
    cout << "-+" << endl;
    cout << "       | ";
    for (int w =0; w < termwdth; w++) cout << " ";
    cout << " |" << endl;
    for (int h = termhght-1; h >= 0; h--) {
        if (h == termhght-1 ) {
            cout << setw(6) << max_val <<  " | ";
        }
        else if (h == termhght/2 ) {
            cout  << setw(6)<< (max_val+min_val)/2. <<  " | ";
        }
        else if (h == 0 ) {
            cout << setw(6) << min_val <<  " | ";
        }
        else {
            cout << "       | ";
        }

        for (int w = termwdth-1; w >= 0; w--) {
            if (matrix [w][h]==1) cout << "@";
            else if (matrix [w][h]==2) cout << ":";
            else cout << " ";
        }

        cout << " |" <<  endl;
    }


    // bottom bar
    cout << "       | ";
    for (int w =0; w < termwdth; w++) cout << " ";
    cout << " |" << endl;

    cout << "       +-";
    for (int w =0; w < termwdth; w++) cout << "-";
    cout << "-+" << endl;
    cout << "                                                                                           " << endl;
    cout << "      CSF                                 cortical depth ->                              WM" << endl;

    cout << endl;
}
    //if (!use_outpath) fout = fin;
   // save_output_nifti(fout, "output1", act, true, use_outpath);
   // save_output_nifti(fout, "output2", layers, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
