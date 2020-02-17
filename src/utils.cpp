
#include "../src/common.h"
#include "../src/utils.h"

// Definitions
void log_welcome(const char* programname) {
    cout << "============="<< endl;
    cout << "LAYNII v1.2.0"<< endl;
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
