
#include "./laynii_lib.h"

// ============================================================================
// Command-line log messages
// ============================================================================

void log_welcome(const char* programname) {
    cout << "======================="<< endl;
    cout << "LayNii v2.5.2          "<< endl;
    cout << "======================="<< endl;
    cout << programname << "\n" << endl;
}

void log_output(const char* filename) {
    cout << "    Writing output as:" << endl;
    cout << "      " << filename << endl;
}

void log_nifti_descriptives(nifti_image* nii) {
    // Print nifti descriptives to command line for debugging
    cout << "    File name: " << nii->fname << endl;
    cout << "    Image details: " << nii->nx  << " X | " << nii->ny << " Y | " << nii->nz << " Z | " << nii->nt << " T " << endl;
    cout << "    Voxel size = " << nii->pixdim[1] << " x " << nii->pixdim[2]   << " x " << nii->pixdim[3] << endl;
    cout << "    Datatype = " << nii->datatype << "\n" << endl;
}

// ============================================================================
// Statistics functions
// ============================================================================

double ren_average(double arr[], int size) {
    int i;
    double sum = 0;
    double avg;

    for (i = 0; i < size; ++i) {
        sum += arr[i];
    }
    avg = double(sum) / size;
    if (avg != avg ) avg = 0; 
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
    if (sqrt( sum) != sqrt( sum) ) sum = 0;
    return sqrt( sum);
}

double ren_correl(double arr1[], double arr2[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double correl = 0;
    double mean1 = ren_average(arr1, size);
    double mean2 = ren_average(arr2, size);

    for (i = 0; i < size; ++i) {
        sum1 += (arr1[i] - mean1) * (arr2[i] - mean2);
        sum2 += (arr1[i] - mean1) * (arr1[i] - mean1);
        sum3 += (arr2[i] - mean2) * (arr2[i] - mean2);
    }
    correl = sum1 / sqrt(sum2 * sum3); 
    if (correl != correl ) correl = 0;
    return correl;
}

double ren_skew(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double mean = ren_average(arr, size);
    double skew = 0; 

    for (i = 0; i < size; ++i) {
        sum1 += (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean);
        sum2 += (arr[i] - mean) * (arr[i] - mean);
    }
    skew = (1 / ((double)size) * sum1) / (pow(1 / ((double)size - 1) * sum2, 1.5)); 
    if (skew != skew ) skew = 0; 
    return (skew);
}

double ren_kurt(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double kurt = 0; 
    double mean = ren_average(arr, size);

    for (i = 0; i < size; ++i) {
        sum1 += (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean) / ((double)size);
        sum2 += (arr[i] - mean) * (arr[i] - mean) / ((double)size);
    }
    kurt =  sum1 / (sum2 * sum2) - 3 ; 
    if (kurt != kurt ) kurt = 0 ;
    return kurt;
}

double ren_autocor(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double autocorr = 0;
    double mean = ren_average(arr, size);

    for (i = 1; i < size; ++i) {
        sum1 += (arr[i]-mean)*(arr[i-1]-mean);
    }
    for (i = 0; i < size; ++i) {
        sum2 += (arr[i]-mean)*(arr[i]-mean);
    }
    autocorr = sum1/sum2 ; 
    if (autocorr != autocorr )  autocorr = 0 ;
    return autocorr;
}


//int ren_add_if_new(int arr[], int size){
//	bool is_new = false :
//	for (i = 1; i < size; ++i) {
//        sum1 += (arr[i]-mean)*(arr[i-1]-mean);
//    }
//}

int ren_most_occurred_number(int nums[], int size){
		int returnvalue = 0;
		int max_count = 0;
		//cout << "\nMost occurred number: ";
		for (int i=0; i<size; i++){
			int count=1;
			for (int j=i+1;j<size;j++)
			if (nums[i]==nums[j])
				count++;
			if (count>max_count)
				max_count = count;
		}

	for (int i=0;i<size;i++){
		int count=1;
			for (int j=i+1;j<size;j++)
				if (nums[i]==nums[j])
					count++;
				if (count==max_count)
					//cout << nums[i] << endl;
					returnvalue =  nums[i]; 
		}
     return returnvalue ; 
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt(pow((x1 - x2) * dX, 2) + pow((y1 - y2) * dY, 2)
                + pow((z1 - z2) * dZ, 2));
}

float dist2d(float x1, float y1, float x2, float y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

float angle(float a, float b, float c) {
    if (a * a + b * b - c * c <= 0) {
        return 3.141592;
    } else {
        return acos((a * a + b * b - c * c) / (2. * a * b));
    }
}

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592))
           * exp(-0.5 * distance * distance / (sigma * sigma));
}

// ============================================================================
// Utility functions
// ============================================================================

void save_output_nifti(const string path, const string tag,  nifti_image* nii,
                       const bool log, const bool use_outpath) {
    ///////////////////////////////////////////////////////////////////////////
    // Note:
    // - 1st argument is the string of the output file name
    //       if there is no explicit output path given, this will be the file
    //       name of the main input data
    //       if there is an explicit output file name given, this will be the
    //       user-defined name following the -output
    //       (including the path and including the file extension)
    // - 2nd argument is the output file name tag, that will be added to the
    //       above argument, this field is ignored, when the flag "use_outpath"
    //       (last argument) is selected.
    // - 3rd argument is the pointer to the data set that is supposed to be
    //       written
    // - 4th argument states if, during the execution of the program an the
    //   writing process should be logged
    //       this argument is optional with the default: TRUE
    // - 5th argument states if the output tag (second argument) should be
    //   ignored or not. This argument is optional the default: FALSE
    //
    // example: save_output_nifti(fout, "VASO_LN", nii_boco_vaso, true, use_outpath);
    ///////////////////////////////////////////////////////////////////////////

    string path_out;

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
            ext = file.substr(pos2);
        } else {  // Determine default extension when no extension given
            basename = file;
            ext = ".nii";
        }

        // Prepare output path
        path_out = dir + sep + basename + "_" + tag + ext;
    }

    // Save nifti
    nifti_set_filenames(nii, path_out.c_str(), 1, 1);
    nifti_image_write(nii);
    if (log) {
        log_output(path_out.c_str());
    }
}


nifti_image* copy_nifti_as_float32_with_scl_slope_and_scl_inter(nifti_image* nii) {
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

    //  Incorporate scaling (scl_slope) and translation (scl_inter) headers
    cout << "  nifti header 'scl slope': " << nii->scl_slope <<endl;
    cout << "  nifti header 'scl inter': " << nii->scl_inter <<endl;
    if (nii->scl_slope != 0) {
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) *= nii->scl_slope;
            *(nii_new_data + i) += nii->scl_inter;
        }
        nii_new->scl_slope = 1.;
        nii_new->scl_inter = 0.;
    } else {
        cout << endl;
        cout << "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "  CAUTION: Nifti scaling parameter 'scl slope' is 0." << endl;
        cout << "    Make sure that your nifti headers are correct!  " << endl;
        cout << "    This program will continue by assuming:         " << endl;
        cout << "      'scl slope = 1' instead of 'scl slope = 0'.   " << endl;
        cout << "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << endl;
    }
    return nii_new;
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

nifti_image* copy_nifti_as_double(nifti_image* nii) {
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
    nii_new->nbyper = sizeof(double);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    double* nii_new_data = static_cast<double*>(nii_new->data);

    // NOTE(Faruk): See nifti1.h for notes on data types
    // ------------------------------------------------------------------------
    if (nii->datatype == 2) {  // NIFTI_TYPE_UINT8
        uint8_t* nii_data = static_cast<uint8_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 512) {  // NIFTI_TYPE_UINT16
        uint16_t* nii_data = static_cast<uint16_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 768) {  // NIFTI_TYPE_UINT32
        uint32_t* nii_data = static_cast<uint32_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 1280) {  // NIFTI_TYPE_UINT64
        uint64_t* nii_data = static_cast<uint64_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 256) {  // NIFTI_TYPE_INT8
        int8_t* nii_data = static_cast<int8_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 4) {  // NIFTI_TYPE_INT16
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 8) {  // NIFTI_TYPE_INT32
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 1024) {  // NIFTI_TYPE_INT64
        int64_t* nii_data = static_cast<int64_t*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 16) {  // NIFTI_TYPE_FLOAT32
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else if (nii->datatype == 64) {  // NIFTI_TYPE_FLOAT64
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nii_new->nvox; ++i) {
            *(nii_new_data + i) = static_cast<double>(*(nii_data + i));
        }
    } else {
        cout << "Warning! Unrecognized nifti data type!" << endl;
    }

    // Replace nans with zeros
    for (int i = 0; i < nii->nvox; ++i) {
        if (*(nii_new_data + i)!= *(nii_new_data + i)) {
            *(nii_new_data + i) = 0;
        }
    }

    return nii_new;
}

nifti_image* copy_nifti_as_int32(nifti_image* nii) {
    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_INT32;
    nii_new->nbyper = sizeof(int32_t);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    int nr_voxels = nii_new->nvox;

    int32_t* nii_new_data = static_cast<int32_t*>(nii_new->data);

    // NOTE(Faruk): See nifti1.h for notes on data types
    // ------------------------------------------------------------------------
    if (nii->datatype == 2) {  // NIFTI_TYPE_UINT8
        uint8_t* nii_data = static_cast<uint8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 512) {  // NIFTI_TYPE_UINT16
        uint16_t* nii_data = static_cast<uint16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 768) {  // NIFTI_TYPE_UINT32
        uint32_t* nii_data = static_cast<uint32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 1280) {  // NIFTI_TYPE_UINT64
        uint64_t* nii_data = static_cast<uint64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 256) {  // NIFTI_TYPE_INT8
        int8_t* nii_data = static_cast<int8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 4) {  // NIFTI_TYPE_INT16
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 8) {  // NIFTI_TYPE_INT32
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 1024) {  // NIFTI_TYPE_INT64
        int64_t* nii_data = static_cast<int64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 16) {  // NIFTI_TYPE_FLOAT32
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 64) {  // NIFTI_TYPE_FLOAT64
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int32_t>(*(nii_data + i));
        }
    } else {
        cout << "Warning! Unrecognized nifti data type!" << endl;
    }

    // Replace nans with zeros
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_new_data + i)!=*(nii_new_data + i)) {
            *(nii_new_data + i) = 0;
        }
    }

    return nii_new;
}

nifti_image* copy_nifti_as_float16(nifti_image* nii) {
    nifti_image* nii_new = nifti_copy_nim_info(nii);
    // NOTE(Renzo): I know that it is not suppoded to look like INT. This is an
    // unlucky naming convention. It is a float16, trust me.
    nii_new->datatype = NIFTI_TYPE_INT16;
    nii_new->nbyper = sizeof(short);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    int nr_voxels = nii_new->nvox;

    short *nii_new_data = static_cast<short*>(nii_new->data);

    // NOTE(Faruk): See nifti1.h for notes on data types
    // ------------------------------------------------------------------------
    if (nii->datatype == 2) {  // NIFTI_TYPE_UINT8
        uint8_t* nii_data = static_cast<uint8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 512) {  // NIFTI_TYPE_UINT16
        uint16_t* nii_data = static_cast<uint16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 768) {  // NIFTI_TYPE_UINT32
        uint32_t* nii_data = static_cast<uint32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 1280) {  // NIFTI_TYPE_UINT64
        uint64_t* nii_data = static_cast<uint64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 256) {  // NIFTI_TYPE_INT8
        int8_t* nii_data = static_cast<int8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 4) {  // NIFTI_TYPE_INT16
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 8) {  // NIFTI_TYPE_INT32
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 1024) {  // NIFTI_TYPE_INT64
        int64_t* nii_data = static_cast<int64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 16) {  // NIFTI_TYPE_FLOAT32
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else if (nii->datatype == 64) {  // NIFTI_TYPE_FLOAT64
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = (short)((double) (*(nii_data + i) * 1000));
        }
    } else {
        cout << "Warning! Unrecognized nifti data type!" << endl;
    }
    nii_new->scl_slope = nii->scl_slope / 1000.;
    // Replace nans with zeros
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nii_new_data + i)!=*(nii_new_data + i)) {
            *(nii_new_data + i) = 0;
        }
    }

    return nii_new;
}

nifti_image* copy_nifti_as_int16(nifti_image* nii) {
    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_INT16;
    nii_new->nbyper = sizeof(int16_t);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    int nr_voxels = nii_new->nvox;

    int16_t* nii_new_data = static_cast<int16_t*>(nii_new->data);

    // NOTE(Faruk): See nifti1.h for notes on data types
    // ------------------------------------------------------------------------
    if (nii->datatype == 2) {  // NIFTI_TYPE_UINT8
        uint8_t* nii_data = static_cast<uint8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 512) {  // NIFTI_TYPE_UINT16
        uint16_t* nii_data = static_cast<uint16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 768) {  // NIFTI_TYPE_UINT32
        uint32_t* nii_data = static_cast<uint32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 1280) {  // NIFTI_TYPE_UINT64
        uint64_t* nii_data = static_cast<uint64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 256) {  // NIFTI_TYPE_INT8
        int8_t* nii_data = static_cast<int8_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 4) {  // NIFTI_TYPE_INT16
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 8) {  // NIFTI_TYPE_INT32
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 1024) {  // NIFTI_TYPE_INT64
        int64_t* nii_data = static_cast<int64_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 16) {  // NIFTI_TYPE_FLOAT32
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
        }
    } else if (nii->datatype == 64) {  // NIFTI_TYPE_FLOAT64
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<int16_t>(*(nii_data + i));
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


// ============================================================================
// Faruk's favorite functions
// ============================================================================
std::tuple<uint32_t, uint32_t, uint32_t> ind2sub_3D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y) {

    uint32_t z = linear_index / (size_x * size_y);
    uint32_t temp = linear_index % (size_x * size_y);
    uint32_t y = temp / size_x;
    uint32_t x = temp % size_x;

    return std::make_tuple(x, y, z);
}

std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> ind2sub_4D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y,
    const uint32_t size_z) {

    uint32_t t = linear_index / (size_x * size_y * size_z);
    uint32_t temp = linear_index % (size_x * size_y * size_z);
    uint32_t z = temp / (size_x * size_y);
    temp = linear_index % (size_x * size_y);
    uint32_t y = temp / size_x;
    uint32_t x = temp % size_x;

    return std::make_tuple(x, y, z, t);
}


uint32_t sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                    const uint32_t size_x, const uint32_t size_y) {
    return size_x * size_y * z + size_x * y + x;
}

uint32_t sub2ind_4D(const uint32_t x, const uint32_t y, const uint32_t z, const uint32_t t,
                    const uint32_t size_x, const uint32_t size_y, const uint32_t size_z) {
    return size_x * size_y * size_z * t + size_x * size_y * z + size_x * y + x;
}


std::tuple<float, float> simplex_closure_2D(float x, float y) {
    float component_sum = x + y;
    float x_new = x / component_sum;
    float y_new = y / component_sum;
    return std::make_tuple(x_new, y_new);
}

std::tuple<float, float> simplex_perturb_2D(float x, float y, float a, float b) {
    float x_new = x * a;
    float y_new = y * b;
    tie(x_new, y_new) = simplex_closure_2D(x_new, y_new);
    return std::make_tuple(x_new, y_new);
}

// std::tuple<float, float> simplex_power_2D(float x, float y, float a) {
//     float x_new = std::pow(x, a);
//     float y_new = std::pow(y, a);
//     tie(x_new, y_new) = simplex_closure_2D(x_new, y_new);
//     return std::make_tuple(x_new, y_new);
// }

// ============================================================================
// Smoothing
// ============================================================================
nifti_image* iterative_smoothing(nifti_image* nii_in, int iter_smooth,
                                 nifti_image* nii_mask, int32_t mask_value) {

    // Copy input niftis TODO[Faruk]: Understand why I have to do this
    nifti_image* temp1 = copy_nifti_as_float32(nii_in);
    nifti_image* temp2 = copy_nifti_as_int32(nii_mask);

    float* nii_in_data = static_cast<float*>(temp1->data);
    int32_t* nii_mask_data = static_cast<int32_t*>(temp2->data);

    // Get dimensions of input
    const uint32_t size_x = temp1->nx;
    const uint32_t size_y = temp1->ny;
    const uint32_t size_z = temp1->nz;
    const uint32_t size_t = temp1->nt;
    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;
    const float dX = temp1->pixdim[1];
    const float dY = temp1->pixdim[2];
    const float dZ = temp1->pixdim[3];

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain voxel visits
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_mask_data + i) != 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_mask_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ------------------------------------------------------------------------
    // Prepare output nifti
    nifti_image* nii_smooth = copy_nifti_as_float32(nii_in);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_smooth_data + i) = 0;
    }

    // Pre-compute weights
    float FWHM_val = 1;  // TODO(Faruk): Might tweak this one
    float w_0 = gaus(0, FWHM_val);
    float w_dX = gaus(dX, FWHM_val);
    float w_dY = gaus(dY, FWHM_val);
    float w_dZ = gaus(dZ, FWHM_val);

    uint32_t ix, iy, iz, j;
    for (uint16_t t = 0; t != size_t; ++t) {  // Over 4th dim (e.g. timepoints)
        for (uint16_t n = 0; n != iter_smooth; ++n) {
            cout << "\r    Iteration: " << n+1 << "/" << iter_smooth << flush;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                uint32_t i = *(voi_id + ii);

                if (*(nii_mask_data + i) == mask_value) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    float new_val = 0, total_weight = 0;

                    // Start with the voxel itself
                    new_val += *(nii_in_data + nr_voxels * t + i) * w_0;
                    total_weight += w_0;

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dX;
                            total_weight += w_dX;
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dX;
                            total_weight += w_dX;
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dY;
                            total_weight += w_dY;
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dY;
                            total_weight += w_dY;
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dZ;
                            total_weight += w_dZ;
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_mask_data + j) == mask_value) {
                            new_val += *(nii_in_data + nr_voxels * t + j) * w_dZ;
                            total_weight += w_dZ;
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    // TODO

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    // TODO

                    *(nii_smooth_data + nr_voxels * t + i) = new_val / total_weight;
                }
            }
            // Swap image data for the next iteration
            for (uint32_t i = 0; i != nr_voxels * size_t; ++i) {
                *(nii_in_data + i) = *(nii_smooth_data + i);
            }
        }
        cout << endl;
    }
    return nii_smooth;
}
