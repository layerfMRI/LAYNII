#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_DOWNSAMPLE: Very simple 2X downsampling. It averages 8 voxels neighboring\n"
    "                voxels into one voxel. The important detail of this program is\n"
    "                that it also adjusts the affine header information according to\n"
    "                the downsamling so that when original and downsampled images are\n"
    "                loaded (e.g. in ITKSNAP) they do not appear to shift.\n"
    "\n"
    "Usage:\n"
    "    LN2_DOWNSAMPLE -input anat.nii.gz\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : A 3D nifti image.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;

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
            fin1 = argv[ac];
            fout = argv[ac];
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
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_DOWNSAMPLE");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint64_t size_x = static_cast<uint64_t>(nii1->nx);
    const uint64_t size_y = static_cast<uint64_t>(nii1->ny);
    const uint64_t size_z = static_cast<uint64_t>(nii1->nz);
    const uint64_t size_time = static_cast<uint64_t>(nii1->nt);
    const uint64_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];

    // Prepare dimensions of output
    const uint64_t out_size_x = size_x/2;
    const uint64_t out_size_y = size_y/2;
    const uint64_t out_size_z = size_z/2;
    const uint64_t out_nr_voxels = out_size_z * out_size_y * out_size_x;

    const float out_dX = dX*2;
    const float out_dY = dY*2;
    const float out_dZ = dZ*2;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocating new nifti for downsampled images
    nifti_image* nii_out = nifti_copy_nim_info(nii1);
    nii_out->datatype = NIFTI_TYPE_INT32;
    // nii_out->dim[0] = 4;  // For proper 4D nifti
    nii_out->dim[1] = out_size_x;
    nii_out->dim[2] = out_size_y;
    nii_out->dim[3] = out_size_z;
    nii_out->dim[4] = size_time;
    nii_out->pixdim[1] = out_dX;
    nii_out->pixdim[2] = out_dY;
    nii_out->pixdim[3] = out_dZ;
    nifti_update_dims_from_array(nii_out);
    nii_out->nvox = out_nr_voxels * size_time;
    nii_out->nbyper = sizeof(int32_t);
    nii_out->data = calloc(nii_out->nvox, nii_out->nbyper);
    nii_out->scl_slope = 1;
    int32_t* nii_out_data = static_cast<int32_t*>(nii_out->data);

    // ------------------------------------------------------------------------
    // Shift affine translation terms by half voxel.
    // ------------------------------------------------------------------------
    // NOTE[Faruk]: This is needed to display the downsampled image without
    // apperaing shifted when loaded as `additional image` in ITKSNAP.
    nii_out->sto_xyz.m[0][3] += dX/2; // affine[0, 3]
    nii_out->sto_xyz.m[1][3] += dY/2; // affine[1, 3]
    nii_out->sto_xyz.m[2][3] += dZ/2; // affine[2, 3]

    for (int i = 0; i != out_nr_voxels*size_time; ++i) {
        *(nii_out_data + i) = 0;
    }

    // ========================================================================
    // Find connected clusters
    // ========================================================================
    cout << "  Downsampling..." << endl;

    for (uint64_t iz = 0; iz != size_z-2; iz += 2) {
        for (uint64_t iy = 0; iy != size_y-2; iy += 2) {
            for (uint64_t ix = 0; ix != size_x-2; ix += 2) {

                // Indices of neighboring voxels
                uint64_t v1 = sub2ind_3D_64(ix  , iy  , iz  , size_x, size_y);
                uint64_t v2 = sub2ind_3D_64(ix+1, iy  , iz  , size_x, size_y);
                uint64_t v3 = sub2ind_3D_64(ix  , iy+1, iz  , size_x, size_y);
                uint64_t v4 = sub2ind_3D_64(ix+1, iy+1, iz  , size_x, size_y);
                uint64_t v5 = sub2ind_3D_64(ix  , iy  , iz+1, size_x, size_y);
                uint64_t v6 = sub2ind_3D_64(ix+1, iy  , iz+1, size_x, size_y);
                uint64_t v7 = sub2ind_3D_64(ix  , iy+1, iz+1, size_x, size_y);
                uint64_t v8 = sub2ind_3D_64(ix+1, iy+1, iz+1, size_x, size_y);

                // Data
                float d1 = *(nii_input_data + v1);
                float d2 = *(nii_input_data + v2);
                float d3 = *(nii_input_data + v3);
                float d4 = *(nii_input_data + v4);
                float d5 = *(nii_input_data + v5);
                float d6 = *(nii_input_data + v6);
                float d7 = *(nii_input_data + v7);
                float d8 = *(nii_input_data + v8);

                // Average
                uint64_t v_new = sub2ind_3D_64(ix/2, iy/2, iz/2, out_size_x, out_size_y);
                *(nii_out_data + v_new) = (d1+d2+d3+d4+d5+d6+d7+d8) / 8;
            }
        }
    }
    cout << endl;

    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "downsampled-2X", nii_out, true);

    cout << "\n  Finished." << endl;
    return 0;
}
