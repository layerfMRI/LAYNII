#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PATCH_UNFLATTEN: Inverse of LN2_PATCH_FLATTEN. This program backwards projects\n"
    "                     flattened image values to the folded image."
    "\n"
    "Usage:\n"
    "    LN2_PATCH_UNFLATTEN -values labels.nii -coord_xyz foldedcoords.nii\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -values    : Nifti image with values that will be backwards projected.\n"
    "                 onto the initial folded image. For example a segmentation\n"
    "                 done on the flat image or a flat domain processed map.\n"
    "    -coord_xyz : A 4D nifti file that contains 3D (XYZ) coordinates.\n"
    "                 For example LN2_PATCH_FLATTEN output named 'foldedcoords'.\n"
    "    -ref       : A nifti image that will be used to extract the folded space\n"
    "                 data dimension information. For instance, '-values' input\n"
    "                 to LN2_PATCH_FLATTEN.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n"
    "Note:\n"
    "    - This program is limited to 3D to 3D projections for now."
    "\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL;
    int ac;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-coord_xyz")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_xyz\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-ref")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -ref\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-values'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-coords_uv'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-ref'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }
    nii3 = nifti_image_read(fin3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }

    log_welcome("LN2_PATCH_UNFLATTEN");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

    // Get dimensions of source input
    const int size_x_flat = nii1->nx;
    const int size_y_flat = nii1->ny;
    const int size_z_flat = nii1->nz;
    const int nr_voxels_flat = size_z_flat * size_y_flat * size_x_flat;

    // Get dimensions of reference (target space) input
    const int size_x_folded = nii3->nx;
    const int size_y_folded = nii3->ny;
    const int size_z_folded = nii3->nz;
    const int nr_voxels_folded = size_z_folded * size_y_folded * size_x_folded;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* flat = copy_nifti_as_float32(nii1);
    float* flat_data = static_cast<float*>(flat->data);
    nifti_image* coords_xyz = copy_nifti_as_float32(nii2);
    float* coords_xyz_data = static_cast<float*>(coords_xyz->data);

    // ========================================================================
    // Prepare outputs
    // ========================================================================
    nifti_image* folded = copy_nifti_as_float32(nii3);
    float* folded_data = static_cast<float*>(folded->data);

    for (int i = 0; i != nr_voxels_folded; ++i) {
        *(folded_data + i) = 0;
    }

    // NOTE[Faruk]: Need to have a counter for many-to-one mapping cases 
    nifti_image* folded_density = copy_nifti_as_float32(folded);
    float* folded_density_data = static_cast<float*>(folded_density->data);

    // ========================================================================
    // Visit each voxel to check their coordinate
    // ========================================================================
    for (int i = 0; i != nr_voxels_flat; ++i) {

        if (*(flat_data + i) != 0) {

            float x = *(coords_xyz_data + nr_voxels_flat*0 + i);
            float y = *(coords_xyz_data + nr_voxels_flat*1 + i);
            float z = *(coords_xyz_data + nr_voxels_flat*2 + i);

            // Cast to integer (floor & cast)
            int cell_idx_x = static_cast<int>(x);
            int cell_idx_y = static_cast<int>(y);
            int cell_idx_z = static_cast<int>(z);

            // Folded image cell index
            int j = size_x_folded * size_y_folded * cell_idx_z + size_x_folded * cell_idx_y + cell_idx_x;

            // Write visited voxel value to folded cell
            *(folded_data + j) += *(flat_data + i);
            *(folded_density_data + j) += 1;
        }
    }

    // Take the mean of each projected cell value
    for (int i = 0; i != nr_voxels_folded; ++i) {
        if (*(folded_density_data + i) > 1) {
            *(folded_data + i) /= *(folded_density_data + i);
        }
    }

    save_output_nifti(fout, "backprojected", folded, true);
    save_output_nifti(fout, "backprojected_density", folded_density, true);

    cout << "\n  Finished." << endl;
    return 0;
}
