#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_DOWNSAMPLE: Very simple integer multiple downsampling. Averages an integer\n"
    "                number of neighboring voxels into one voxel."
    "\n"
    "Usage:\n"
    "    LN2_DOWNSAMPLE -input anat.nii.gz -factor 2\n"
    "\n"
    "Options:\n"
    "    -help    : Show this help.\n"
    "    -input   : A 3D nifti image.\n"
    "    -factor  : Downsampling factor. Only positive integers are allowed.\n"
    "    -inplace : (Optional) Save downsampled image in the original resolution.\n"
    "    -output  : (Optional) Output basename for all outputs.\n"
    "    -debug   : (Optional) Save extra intermediate outputs.\n"
    "\n"
    "Notes:\n"
    "    - The important detail of this program is that it also adjusts the affine\n"
    "      header information according to the downsamling so that when original and\n"
    "      downsampled images are loaded (e.g. in ITKSNAP) they do not appear to shift.\n"
    "    - Header is not changed in '-inplace' mode.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    uint64_t factor = 2; 
    bool mode_inplace = false;

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
        } else if (!strcmp(argv[ac], "-factor")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -factor\n");
            } else {
                factor = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-inplace")) {
            mode_inplace = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Check input validity
    if (factor < 1) {
        fprintf(stderr, "  Downsampling factor is not a positive integer.'\n");
        return 2;
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
    const uint64_t out_size_x = size_x/factor;
    const uint64_t out_size_y = size_y/factor;
    const uint64_t out_size_z = size_z/factor;
    const uint64_t out_nr_voxels = out_size_z * out_size_y * out_size_x;

    const float out_dX = dX*static_cast<float>(factor);
    const float out_dY = dY*static_cast<float>(factor);
    const float out_dZ = dZ*static_cast<float>(factor);

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Prepare output
    nifti_image* nii_out;
    float* nii_out_data;

    if (mode_inplace) {
        nii_out = copy_nifti_as_float32(nii_input);
        nii_out_data = static_cast<float*>(nii_out->data);

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_out_data + i) = 0;
        }

    } else {
        // Allocating new nifti for downsampled images
        nii_out = nifti_copy_nim_info(nii_input);
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
        nii_out->nbyper = sizeof(float);
        nii_out->data = calloc(nii_out->nvox, nii_out->nbyper);
        nii_out->scl_slope = 1;
        nii_out_data = static_cast<float*>(nii_out->data);        

        // ------------------------------------------------------------------------
        // Shift affine translation terms by half voxel.
        // ------------------------------------------------------------------------
        // NOTE[Faruk]: This is needed to display the downsampled image without
        // apperaing shifted when loaded as `additional image` in ITKSNAP.
        nii_out->sto_xyz.m[0][3] += dX/static_cast<float>(factor); // affine[0, 3]
        nii_out->sto_xyz.m[1][3] += dY/static_cast<float>(factor); // affine[1, 3]
        nii_out->sto_xyz.m[2][3] += dZ/static_cast<float>(factor); // affine[2, 3]
    }

    // ========================================================================
    // Find connected clusters
    // ========================================================================
    cout << "  Downsampling with a factor of " << factor << "..." << endl;

    for (uint64_t iz = 0; iz < size_z; iz += factor) {
        for (uint64_t iy = 0; iy < size_y; iy += factor) {
            for (uint64_t ix = 0; ix < size_x; ix += factor) {

                float d_new = 0;

                // Collect data from neighboring voxels
                for (uint64_t jz = 0; jz != factor; ++jz) {
                    for (uint64_t jy = 0; jy != factor; ++jy) {
                        for (uint64_t jx = 0; jx != factor; ++jx) {
                            uint64_t v_temp = sub2ind_3D_64(ix+jx, iy+jy, iz+jz, size_x, size_y);
                            d_new += *(nii_input_data + v_temp);
                        }
                    }
                }

                // Average
                d_new /= static_cast<float>(factor*factor*factor);

                if (mode_inplace) {
                // Write new data back to neighboring voxels
                    for (uint64_t jz = 0; jz != factor; ++jz) {
                        for (uint64_t jy = 0; jy != factor; ++jy) {
                            for (uint64_t jx = 0; jx != factor; ++jx) {
                                uint64_t v_temp = sub2ind_3D_64(ix+jx, iy+jy, iz+jz, size_x, size_y);
                                *(nii_out_data + v_temp) = d_new;
                            }
                        }
                    }
                } else {
                    uint64_t v_new = sub2ind_3D_64(ix/factor, iy/factor, iz/factor, out_size_x, out_size_y);
                    *(nii_out_data + v_new) = d_new;
                }
            }
        }
    }
    cout << endl;

    cout << "  Saving output..." << endl;
    std::ostringstream tag;
    tag << factor;
    if (mode_inplace) {
        save_output_nifti(fout, "downsampled_inplace-" + tag.str() + "X", nii_out, true);        
    } else {
        save_output_nifti(fout, "downsampled-" + tag.str() + "X", nii_out, true);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
