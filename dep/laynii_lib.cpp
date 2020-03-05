
#include "./laynii_lib.h"

// ============================================================================
// ============================================================================
// ============================================================================
// Statistics functions

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

double ren_correl(double arr1[], double arr2[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double mean1 = ren_average(arr1, size);
    double mean2 = ren_average(arr2, size);

    for (i = 0; i < size; ++i) {
        sum1 += (arr1[i] - mean1) * (arr2[i] - mean2);
        sum2 += (arr1[i] - mean1) * (arr1[i] - mean1);
        sum3 += (arr2[i] - mean2) * (arr2[i] - mean2);
    }
    return sum1 / sqrt(sum2 * sum3);
}

double ren_skew(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double mean = ren_average(arr, size);

    for (i = 0; i < size; ++i) {
        sum1 += (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean);
        sum2 += (arr[i] - mean) * (arr[i] - mean);
    }
    return ((1 / ((double)size) * sum1) / (pow(1 / ((double)size - 1) * sum2, 1.5)));
}

double ren_kurt(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double mean = ren_average(arr, size);

    for (i = 0; i < size; ++i) {
        sum1 += (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean) * (arr[i] - mean) / ((double)size);
        sum2 += (arr[i] - mean) * (arr[i] - mean) / ((double)size);
    }
    return sum1 / (sum2 * sum2) - 3;
}

double ren_autocor(double arr[], int size) {
    int i;
    double sum1 = 0;
    double sum2 = 0;
    double mean = ren_average(arr, size);

    for (i = 1; i < size; ++i) {
        sum1 += (arr[i]-mean)*(arr[i-1]-mean);
    }
    for (i = 0; i < size; ++i) {
        sum2 += (arr[i]-mean)*(arr[i]-mean);
    }
    return sum1/sum2;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt(pow((x1 - x2) * dX, 2) + pow((y1 - y2) * dY, 2)
                + pow((z1 - z2) * dZ, 2));
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
// ============================================================================
// ============================================================================
// Command-line log messages

void log_welcome(const char* programname) {
    cout << "============="<< endl;
    cout << "LAYNII v1.1.0"<< endl;
    cout << "============="<< endl;
    cout << programname << "\n" << endl;
}

void log_output(const char* filename) {
    cout << "  Writing output as:" << endl;
    cout << "    " << filename << endl;
}

void log_nifti_descriptives(nifti_image* nii) {
    // Print nifti descriptives to command line for debugging
    cout << "  File name: " << nii->fname << endl;
    cout << "    Image details: " << nii->nz << " Slices | " << nii->nx
         << " Phase steps | " << nii->ny << " Read steps | " << nii->nt
         << " Time steps " << endl;
    cout << "    Voxel size = " << nii->pixdim[1] << " x " << nii->pixdim[2]
         << " x " << nii->pixdim[3] << endl;
    cout << "    Datatype = " << nii->datatype << "\n" << endl;
}

// ============================================================================
// ============================================================================
// ============================================================================
// Utility functions

void save_output_nifti(string filename, string prefix, nifti_image* nii,
                       bool log) {
    string outfilename = prefix + "_" + filename;
    nifti_set_filenames(nii, outfilename.c_str(), 1, 1);
    nifti_image_write(nii);
    if (log) {
        log_output(outfilename.c_str());
    }
}

nifti_image* copy_nifti_header_as_float(nifti_image* nii) {
    ///////////////////////////////////////////////////////////////////////////
    // NOTE(Renzo): Fixing potential problems with different input datatypes //
    // here, I am loading them in their native datatype and translate them   //
    // to the datatime I like best.                                          //
    ///////////////////////////////////////////////////////////////////////////

    // NOTE(for future reference): Rick's comments:
    // nifti_copy_nim_info(). It will return with data == NULL.
    // If you need the data allocated, memory use would not change once you do
    // so. There is also nifti_make_new_nim()
    // nifti_image* nifti_make_new_nim(const int64_t dims[], int datatype, int data_fill)

    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    return nii_new;
}

nifti_image* copy_nifti_header_as_int(nifti_image* nii) {
    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_INT32;
    nii_new->nbyper = sizeof(int32_t);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    return nii_new;
}

// ============================================================================
// Faruk's favorite functions
std::tuple<uint32_t, uint32_t, uint32_t> ind2sub_3D(
    const uint32_t linear_index, const uint32_t size_x, const uint32_t size_y) {
    uint32_t z = linear_index / (size_x * size_y);
    uint32_t temp = linear_index % (size_x * size_y);
    uint32_t y = temp / size_x;
    uint32_t x = temp % size_x;
    return std::make_tuple(x, y, z);
}

uint32_t sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                    const uint32_t size_x, const uint32_t size_y) {
    return size_x * size_y * z + size_x * y + x;
}
