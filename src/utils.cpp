
#include "../src/common.h"
#include "../src/utils.h"

// Definitions
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

nifti_image* recreate_nii_with_float_datatype(nifti_image* nii) {
    ///////////////////////////////////////////////////////////////////////////
    // NOTE(Renzo): Fixing potential problems with different input datatypes //
    // here, I am loading them in their native datatype and translate them   //
    // to the datatime I like best.                                          //
    ///////////////////////////////////////////////////////////////////////////

    // Get dimensions of input
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int size_z = nii->nz;
    const int size_t = nii->nt;
    const int nx = nii->nx;
    const int nxy = nii->nx * nii->ny;
    const int nxyz = nii->nx * nii->ny * nii->nz;

    // NOTE(for future reference): Rick's comments:
    // nifti_copy_nim_info(). It will return with data == NULL.
    // If you need the data allocated, memory use would not change once you do so.
    // There is also nifti_make_new_nim()
    // nifti_image * nifti_make_new_nim(const int64_t dims[], int datatype, int data_fill)

    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    if (nii->datatype == NIFTI_TYPE_INT8 || nii->datatype == DT_UINT8) {
        int8_t* temp = static_cast<int8_t*>(nii->data);
        FOR_EACH_VOXEL_TZYX
            *(nii_new_data + VOXEL_ID) = static_cast<float>(*(temp + VOXEL_ID));
        END_FOR_EACH_VOXEL_TZYX
    } else if (nii->datatype == NIFTI_TYPE_INT16) {
        int16_t* temp = static_cast<int16_t*>(nii->data);
        FOR_EACH_VOXEL_TZYX
            *(nii_new_data + VOXEL_ID) = static_cast<float>(*(temp + VOXEL_ID));
        END_FOR_EACH_VOXEL_TZYX
    } else if (nii->datatype == NIFTI_TYPE_INT32) {
        int* temp = static_cast<int*>(nii->data);
        FOR_EACH_VOXEL_TZYX
            *(nii_new_data + VOXEL_ID) = static_cast<float>(*(temp + VOXEL_ID));
        END_FOR_EACH_VOXEL_TZYX
    } else if (nii->datatype == NIFTI_TYPE_FLOAT32) {
        float* temp = static_cast<float*>(nii->data);
        FOR_EACH_VOXEL_TZYX
            *(nii_new_data + VOXEL_ID) = static_cast<float>(*(temp + VOXEL_ID));
        END_FOR_EACH_VOXEL_TZYX
    }
    return nii_new;
}
