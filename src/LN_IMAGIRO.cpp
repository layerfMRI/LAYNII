

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_IMAGIRO : Generates a 3D matrix of columns and layers.\n"
    "             It Does the ORIGAMI backwards; It unfolds stuff.\n"
    "\n"
    "Usage:\n"
    "    LN_IMAGIRO -layers equi_dist_layers.nii -columns columnar_coordinated.nii -data data2unfold.nii\n"
    "\n"
    "  Test application in the test_data folder: \n"
    "    ../LN_IMAGIRO -layers sc_layers_3dcolumns.nii -column_file sc_columns_3dcolumns.nii -data sc_BOLD_act.nii \n"
    "\n"
    "  a potiential application case is described on this blog post: \n"
    "    https://layerfmri.com/columns/ \n"
    "    \n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -layers     : Nifti (.nii) file that contains layer.\n"
    "    -columns    : Nifti (.nii) file that contains columns.\n"
    "    -data       : Data that will be unfolded.\n"
    "\n"
    "Notes:\n"
    "    - All inputs should have the same dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char* fin_layers = NULL, * fin_columns = NULL, * fin_data = NULL;
    int ac;
    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -layers\n");
                return 1;
            }
            fin_layers = argv[ac];
        } else if (!strcmp(argv[ac], "-columns")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -columns\n");
                return 1;
            }
            fin_columns = argv[ac];
        } else if (!strcmp(argv[ac], "-data")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -data\n");
                return 1;
            }
            fin_data = argv[ac];
        } else {
            fprintf(stderr, " ** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_columns) {
        fprintf(stderr, " ** missing option '-columns'\n");
        return 1;
    }
    if (!fin_layers) {
        fprintf(stderr, " ** missing option '-layers'\n");
        return 1;
    }
    if (!fin_data) {
        fprintf(stderr, " ** missing option '-data'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nim_column_r = nifti_image_read(fin_columns, 1);
    if (!nim_column_r) {
        fprintf(stderr, " ** failed to read NIfTI from '%s'\n", fin_columns);
        return 2;
    }
    nifti_image* nim_layers_r = nifti_image_read(fin_layers, 1);
    if (!nim_layers_r) {
        fprintf(stderr, " ** failed to read NIfTI from '%s'\n", fin_layers);
        return 2;
    }
    nifti_image* nim_data_r = nifti_image_read(fin_data, 1);
    if (!nim_data_r) {
        fprintf(stderr, " ** failed to read NIfTI from '%s'\n", fin_data);
        return 2;
    }

    log_welcome("LN_IMAGIRO");
    log_nifti_descriptives(nim_layers_r);
    log_nifti_descriptives(nim_column_r);
    log_nifti_descriptives(nim_data_r);

    // Get dimensions of input
    int size_z = nim_layers_r->nz;
    int size_x = nim_layers_r->nx;
    int size_y = nim_layers_r->ny;
    int size_time = nim_data_r->nt;
    int nx = nim_layers_r->nx;
    int nxy = nim_layers_r->nx * nim_layers_r->ny;
    int nxyz = nim_layers_r->nx * nim_layers_r->ny * nim_layers_r->nz;
    int nr_voxels = nim_layers_r->nvox;
    float dX = nim_layers_r->pixdim[1];
    float dY = nim_layers_r->pixdim[2];
    float dZ = nim_layers_r->pixdim[3];

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nim_layers = copy_nifti_as_int32(nim_layers_r);
    int32_t* nim_layers_data = static_cast<int32_t*>(nim_layers->data);

    nifti_image* nim_columns = copy_nifti_as_int32(nim_column_r);
    int32_t* nim_columns_data = static_cast<int32_t*>(nim_columns->data);

    nifti_image* nim_data = copy_nifti_as_float32(nim_data_r);
    float* nim_data_data = static_cast<float*>(nim_data->data);

    // ========================================================================
    // Finding number of layers/columns
    int nr_layers = 0, nr_columns = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_layers_data + i) > nr_layers) {
            nr_layers = *(nim_layers_data + i);
        }
        if (*(nim_columns_data + i) > nr_columns) {
            nr_columns = *(nim_columns_data + i);
        }
    }
    cout << "  There are " << nr_layers << " layers " << endl;
    cout << "  There are " << nr_columns << " columns " << endl;

    // ========================================================================
    ////////////////////////////////
    // Allocating necessary files //
    ////////////////////////////////
    nifti_image* imagiro = nifti_copy_nim_info(nim_data);
    imagiro->datatype = NIFTI_TYPE_FLOAT32;
    imagiro->nbyper = sizeof(float);

    imagiro->pixdim[1] = 1;
    imagiro->pixdim[2] = 1;
    imagiro->pixdim[3] = 1;

    imagiro->nx = nr_columns;
    imagiro->ny = size_z;
    imagiro->nz = nr_layers;
    imagiro->nt = size_time;

    imagiro->dim[1] = imagiro->nx;
    imagiro->dim[2] = imagiro->ny;
    imagiro->dim[3] = imagiro->nz;
    imagiro->nvox = imagiro->nx * imagiro->ny * imagiro->nz * imagiro->nt
                    * imagiro->nu * imagiro->nv * imagiro->nw;

    int size_x_imagiro = imagiro->nx;
    int size_y_imagiro = imagiro->ny;
    int size_z_imagiro = imagiro->nz;
    int size_t_imagiro = imagiro->nt;
    int nx_imagiro = imagiro->nx;
    int nxy_imagiro = imagiro->nx * imagiro->ny;
    int nxyz_imagiro = imagiro->nx * imagiro->ny * imagiro->nz;
    float dX_imagiro = imagiro->pixdim[1];
    float dY_imagiro = imagiro->pixdim[2];
    float dZ_imagiro = imagiro->pixdim[3];

    cout << "  imagiro layer dim a " << imagiro->nz << endl;
    // nifti_update_dims_from_array(imagiro) ;
    cout << "  imagiro layer dim b " << imagiro->nz << endl;
    // cout << imagiro->pixdim[1] << endl;
    // cout << imagiro->pixdim[2] << endl;
    // cout << imagiro->pixdim[3] << endl;
    cout << "  imagiro column dim " << imagiro->nx << endl;
    cout << "  imagiro depth dim " << imagiro->ny << endl;
    cout << "  imagiro layer dim " << imagiro->nz << endl;
    cout << "  imagiro time dim " << imagiro->nt << endl;

    imagiro->data = calloc(imagiro->nvox, imagiro->nbyper);
    float* imagiro_data = static_cast<float*>(imagiro->data);

    int nii_ok = 4;
    nifti_nim_has_valid_dims(imagiro, nii_ok);
    nifti_nim_is_valid(imagiro, nii_ok);

    nifti_image* imagiro_vnr = nifti_copy_nim_info(imagiro);
    imagiro_vnr->datatype = NIFTI_TYPE_INT32;
    imagiro_vnr->nbyper = sizeof(int);
    imagiro_vnr->data = calloc(imagiro_vnr->nvox, imagiro_vnr->nbyper);
    int* imagiro_vnr_data = static_cast<int*>(imagiro_vnr->data);

    // ========================================================================
    ////////////////////////////////////////
    // Resample across Layers and columns //
    ////////////////////////////////////////
    cout << "  Resampling..." << endl;
    float value_ofinput_data = 0;
    imagiro->nz = nr_layers;
    imagiro->nx = nr_columns;
    imagiro->ny = size_z;

    for (int iz = 0; iz < size_z_imagiro; ++iz) {
        for (int iy = 0; iy < size_y_imagiro; ++iy) {
            for (int ix = 0; ix < size_x_imagiro; ++ix) {
                *(imagiro_data + nxy_imagiro * iz + nx_imagiro * iy + ix) = 0;
                *(imagiro_vnr_data + nxy_imagiro * iz + nx_imagiro * iy + ix) = 0;
            }
        }
    }

    ///////////////////////////////////////////////////////
    // Calculating the number of voxels per layer/column //
    ///////////////////////////////////////////////////////
    cout << "  Calculating the number of voxels per layer column..." << endl;
    int lay = 0, col = 0, dep = 0;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                if (*(nim_layers_data + voxel_i) > 0
                    && *(nim_columns_data + voxel_i) > 0) {
                    lay = *(nim_layers_data + voxel_i) - 1;
                    col = *(nim_columns_data + voxel_i) - 1;
                    dep = iz;
                    *(imagiro_vnr_data + nxy_imagiro * lay + nx_imagiro * dep + col) += 1;
                }
            }
        }
    }

    //////////////////////////////////////////
    // Averaging all voxels in layer\column //
    //////////////////////////////////////////
    cout << "  Averaging all voxels in layer column..." << endl;
    for (int it = 0; it < size_time; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;

                    if (*(nim_layers_data + voxel_i) > 0
                        && *(nim_columns_data + voxel_i) > 0) {
                        lay = *(nim_layers_data + voxel_i) - 1;
                        col = *(nim_columns_data + voxel_i) - 1;
                        dep = iz;
                        int voxel_j = nxy_imagiro * lay + nx_imagiro * dep + col;

                        value_ofinput_data =
                            *(nim_data_data + nxyz * it + voxel_i)
                            / *(imagiro_vnr_data + voxel_j);
                        *(imagiro_data + nxyz_imagiro * it + voxel_j) += value_ofinput_data;
                    }
                }
            }
        }
    }
    //////////////////
    // Fixing holes //
    //////////////////
    cout << "  Fixing holes..." << size_time << endl;
    int vinc = 2;  // vicinity
    float value_to_fill = 0;  // Average value in vicinity
    float nr_vic = 0;  // Number of non-zero voxels in vicinity

    for (int iz = 0; iz < size_z_imagiro; ++iz) {
        for (int iy = 0; iy < size_y_imagiro; ++iy) {
            for (int ix = 0; ix < size_x_imagiro; ++ix) {
                int voxel_i = nxy_imagiro * iz + nx_imagiro * iy + ix;

                if (*(imagiro_vnr_data + voxel_i) == 0) {
                    for (int it = 0; it < size_time; ++it) {
                        value_to_fill = 0;
                        nr_vic = 0;

                        int jy_start = max(0, iy - vinc);
                        int jy_stop = min(iy + vinc + 1, size_x_imagiro);
                        int jx_start = max(0, ix - vinc);
                        int jx_stop = min(ix + vinc + 1, size_y_imagiro);
                        int jz_start = max(0, iz - vinc);
                        int jz_stop = min(iz + vinc + 1, size_z_imagiro);

                        for (int jz = jz_start; jz < jz_stop; ++jz) {
                            for (int jy = jy_start; jy < jy_stop; ++jy) {
                                for (int jx = jx_start; jx < jx_stop; ++jx) {
                                    int voxel_j = nxy_imagiro * jz + nx_imagiro * jy + jx;

                                    if (*(imagiro_vnr_data + voxel_j) != 0) {
                                        nr_vic += 1;
                                        value_to_fill += static_cast<float>(*(imagiro_data + nxyz_imagiro * it + voxel_j));
                                    }
                                }
                            }
                        }
                        *(imagiro_data + nxyz_imagiro * it + voxel_i) = value_to_fill / nr_vic;
                    }
                }
            }
        }
    }

    ////////////////////////////////////////
    // Taking care of scale factor if any //
    ////////////////////////////////////////
    imagiro->scl_slope = nim_data_r->scl_slope;
    imagiro->iname_offset = nim_data_r->iname_offset;
    if (nim_data_r->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    save_output_nifti(fin_data, "unfolded", imagiro, true);
    save_output_nifti(fin_data, "nr_voxels", imagiro_vnr, true);

    cout << "  Finished." << endl;
    return 0;
}
