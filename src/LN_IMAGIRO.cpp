
#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_IMAGIRO : This program does the opposite of ORIGAMI. It unfolds stuff.\n"
    "\n"
    "    Generates a 3D matrix of columns and layers.\n"
    "\n"
    "Usage:\n"
    "    LN_IMAGIRO -layer_file equi_dist_layers.nii -column_file columnar_coordinated.nii -data data2unfold.nii\n"
    "\n"
    "Options:\n"
    "   -help         : Show this help.\n"
    "   -layer_file   : Nifti (.nii) file that contains layer.\n"
    "   -columns_file : Nifti (.nii) file that contains columns.\n"
    "   -data         : Data that will be unfolded.\n"
    "\n"
    "Notes:\n"
    "    - All inputs should have the same dimensions.\n"
    "    - This program now supports INT16, INT32 and FLOAT32.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * layer_filename = NULL, * fout = NULL, * column_filename = NULL, * data_filename = NULL;
    int ac, twodim = 0, do_masking = 0;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        }
        else if (!strcmp(argv[ac], "-layer_file")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -layer_file\n");
                return 1;
            }
            layer_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-column_file")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -column_file\n");
                return 1;
            }
            column_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-data")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -data\n");
                return 1;
            }
            data_filename = argv[ac];  // Assign pointer, no string copy
        } else {
            fprintf(stderr, " * * invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!column_filename) {
        fprintf(stderr, " * * missing option '-column_file'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_column_r = nifti_image_read(column_filename, 1);
    if (!nim_column_r) {
        fprintf(stderr, " * * failed to read layer NIfTI image from '%s'\n", column_filename);
        return 2;
    }
    if (!layer_filename) {
        fprintf(stderr, " * * missing option '-layer_file'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_layers_r = nifti_image_read(layer_filename, 1);
    if (!nim_layers_r) {
        fprintf(stderr, " * * failed to read layer NIfTI image from '%s'\n", layer_filename);
        return 2;
    }
    if (!data_filename) {
        fprintf(stderr, " * * missing option '-data'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_data_r = nifti_image_read(data_filename, 1);
    if (!nim_data_r) {
        fprintf(stderr, " * * failed to read layer NIfTI image from '%s'\n", data_filename);
        return 2;
    }

    // Get dimensions of input
    int sizeSlice = nim_layers_r->nz;
    int sizePhase = nim_layers_r->nx;
    int sizeRead = nim_layers_r->ny;
    int nx = nim_layers_r->nx;
    int nxy = nim_layers_r->nx * nim_layers_r->ny;
    int nxyz = nim_layers_r->nx * nim_layers_r->ny * nim_layers_r->nz;
    float dX = nim_layers_r->pixdim[1];
    float dY = nim_layers_r->pixdim[2];
    float dZ = nim_layers_r->pixdim[3];
    int nrep = nim_data_r->nt;

    if (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image * nim_layers = nifti_copy_nim_info(nim_layers_r);
    nim_layers->datatype = NIFTI_TYPE_INT16;
    nim_layers->nbyper = sizeof(short);
    nim_layers->data = calloc(nim_layers->nvox, nim_layers->nbyper);
    short * nim_layers_data = (short * ) nim_layers->data;

    nifti_image * nim_columns = nifti_copy_nim_info(nim_column_r);
    nim_columns->datatype = NIFTI_TYPE_INT16;
    nim_columns->nbyper = sizeof(short);
    nim_columns->data = calloc(nim_columns->nvox, nim_columns->nbyper);
    short * nim_columns_data = (short * ) nim_columns->data;
    nim_columns->scl_slope = nim_column_r->scl_slope;

    nifti_image * nim_data = nifti_copy_nim_info(nim_data_r);
    nim_data->datatype = NIFTI_TYPE_FLOAT32;
    nim_data->nbyper = sizeof(float);
    nim_data->data = calloc(nim_data->nvox, nim_data->nbyper);
    float * nim_data_data = (float * ) nim_data->data;

    nim_data->scl_slope = nim_data_r->scl_slope;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_column_r->datatype == NIFTI_TYPE_FLOAT32) {
        float * nim_column_r_data = (float * ) nim_column_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_columns_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_column_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_column_r->datatype == NIFTI_TYPE_INT16) {
        short * nim_column_r_data = (short * ) nim_column_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_columns_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_column_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    // if (* (nim_columns_data + nxyz * it + nxy * islice + nx * ix + iy) > 0) cout << * (nim_columns_data + nxyz * it + nxy * islice + nx * ix + iy) << endl;
                }
            }
        }
    }
    if (nim_column_r->datatype == NIFTI_TYPE_INT32) {
        int * nim_column_r_data = (int * ) nim_column_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_columns_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_column_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_FLOAT32) {
        float * nim_layers_r_data = (float * ) nim_layers_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_INT16) {
        short * nim_layers_r_data = (short * ) nim_layers_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_layers_r->datatype == NIFTI_TYPE_INT32) {
        int * nim_layers_r_data = (int * ) nim_layers_r->data;
        int it = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_layers_data + nxyz * it + nxy * islice + nx * ix + iy) = (short) (*(nim_layers_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_data_r->datatype == NIFTI_TYPE_FLOAT32) {
        cout << "Input data is float " << endl;
        float * nim_data_r_data = (float * ) nim_data_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_data_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_data_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_data_r->datatype == NIFTI_TYPE_INT16) {
        cout << "Input data is int 16 " << endl;
        short * nim_data_r_data = (short * ) nim_data_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        // TODO: ask Rick about this if (* (nim_data_r_data + nxyz * it + nxy * islice + nx * ix + iy) > 0) cout << * (nim_data_r_data + nxyz * it + nxy * islice + nx * ix + iy) << endl;
                        *(nim_data_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_data_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_data_r->datatype == NIFTI_TYPE_INT32) {
        cout << "Input data is int 32 " << endl;
        int * nim_data_r_data = (int * ) nim_data_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_data_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_data_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    cout << " " << sizeSlice << " Slices | " << sizePhase << " Phase_steps | " << sizeRead << " Read_steps | " << nrep << " Time_steps " << endl;
    cout << "  Voxel size = " << dX << " x " << dY << " x " << dZ << endl;
    cout << "  Datatype of Layers mask = " << nim_layers->datatype << endl;
    cout << "  Datatype of Columns mask = " << nim_columns->datatype << endl;

    //////////////////////////////
    // Finding number of layers //
    //////////////////////////////
    int layernumber = 0;
    int columnnumber = 0;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > layernumber) layernumber = *(nim_layers_data + nxy * iz + nx * ix + iy);
                if (*(nim_columns_data + nxy * iz + nx * ix + iy) > columnnumber) {
                    columnnumber = *(nim_columns_data + nxy * iz + nx * ix + iy);
                    //cout << " iz " << iz << " ix " << ix << " iy " << iy << " " << columnnumber << " " << nim_columns->scl_slope << endl;
                }
            }
        }
    }
    cout << "  There are " << layernumber << " layers " << endl;
    cout << "  There are " << columnnumber << " columns " << endl;

    // cout << "  nim_data_data->cal_max " << * nifti_image_to_ascii(nim_data) << endl;

    ////////////////////////////////
    // Allocating necessary files //
    ////////////////////////////////
    nifti_image * imagiro = nifti_copy_nim_info(nim_data);
    imagiro->datatype = NIFTI_TYPE_FLOAT32;
    imagiro->nbyper = sizeof(float);

    imagiro->pixdim[1] = 1;
    imagiro->pixdim[2] = 1;
    imagiro->pixdim[3] = 1;

    imagiro->nx = columnnumber;
    imagiro->ny = sizeSlice;
    imagiro->nz = layernumber;
    imagiro->nt = nrep;

    imagiro->dim[1] = imagiro->nx;
    imagiro->dim[2] = imagiro->ny;
    imagiro->dim[3] = imagiro->nz;
    imagiro->nvox = imagiro->nx * imagiro->ny * imagiro->nz * imagiro->nt * imagiro->nu * imagiro->nv * imagiro->nw;

    int sizeColumn_imagiro = imagiro->nx;
    int sizeDepth_imagiro = imagiro->ny;
    int sizeLayer_imagiro = imagiro->nz;
    int nrep_imagiro = imagiro->nt;
    int nx_imagiro = imagiro->nx;
    int nxy_imagiro = imagiro->nx * imagiro->ny;
    int nxyz_imagiro = imagiro->nx * imagiro->ny * imagiro->nz;
    float dX_imagiro = imagiro->pixdim[1];
    float dY_imagiro = imagiro->pixdim[2];
    float dZ_imagiro = imagiro->pixdim[3];

    cout << "  imagiro layer dim a " << imagiro->nz << endl;
    //nifti_update_dims_from_array(imagiro) ;
    cout << "  imagiro layer dim b " << imagiro->nz << endl;
    //  cout << imagiro->pixdim[1] << endl;
    //  cout << imagiro->pixdim[2] << endl;
    //  cout << imagiro->pixdim[3] << endl;
    cout << "  imagiro column dim " << imagiro->nx << endl;
    cout << "  imagiro depth dim " << imagiro->ny << endl;
    cout << "  imagiro layer dim " << imagiro->nz << endl;
    cout << "  imagiro time dim " << imagiro->nt << endl;

    imagiro->data = calloc(imagiro->nvox, imagiro->nbyper);
    float * imagiro_data = (float * ) imagiro->data;

    int nii_ok = 4;
    nifti_nim_has_valid_dims(imagiro, nii_ok);
    nifti_nim_is_valid(imagiro, nii_ok);

    nifti_image * imagiro_vnr = nifti_copy_nim_info(imagiro);
    imagiro_vnr->datatype = NIFTI_TYPE_INT32;
    imagiro_vnr->nbyper = sizeof(int);
    imagiro_vnr->data = calloc(imagiro_vnr->nvox, imagiro_vnr->nbyper);
    int * imagiro_vnr_data = (int * ) imagiro_vnr->data;

    // cout << "  nii ok = " << nii_ok << endl;

    ////////////////////////////////////////
    // Resample across Layers and columns //
    ////////////////////////////////////////
    cout << "resampling " << endl;
    short layer_zi = 0;
    short column_yi = 0;
    short depth_xi = 0;
    float value_ofinput_data = 0;
    imagiro->nz = layernumber;
    imagiro->nx = columnnumber;
    imagiro->ny = sizeSlice;

    for (int iz = 0; iz < sizeLayer_imagiro; ++iz) {
        for (int iy = 0; iy < sizeColumn_imagiro; ++iy) {
            for (int ix = 0; ix < sizeDepth_imagiro; ++ix) {
                // cout << iz << " " << ix << " " << iy << endl;
                *(imagiro_data + nxy_imagiro * iz + nx_imagiro * ix + iy) = 0;
                *(imagiro_vnr_data + nxy_imagiro * iz + nx_imagiro * ix + iy) = 0;
                // *(nim_columns_data + nxy * iz + nx * ix + iy) = 1 ;
            }
        }
    }

    ///////////////////////////////////////////////////////
    // Calculating the number of voxels per layer/column //
    ///////////////////////////////////////////////////////

    cout << "  Calculating the number of voxels per layer column..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_layers_data + nxy * iz + nx * ix + iy) > 0 && *(nim_columns_data + nxy * iz + nx * ix + iy) > 0) {
                    layer_zi = *(nim_layers_data + nxy * iz + nx * ix + iy)-1;
                    column_yi = *(nim_columns_data + nxy * iz + nx * ix + iy)-1;
                    depth_xi = iz;
                    *(imagiro_vnr_data + nxy_imagiro * layer_zi + nx_imagiro * depth_xi + column_yi) = *(imagiro_vnr_data + nxy_imagiro * layer_zi + nx_imagiro * depth_xi + column_yi)+ 1;
                }
            }
        }
    }

    //////////////////////////////////////////
    // Averaging all voxels in layer\column //
    //////////////////////////////////////////
    cout << "  Averaging all voxels in layer column..." << endl;
    for (int it = 0; it < nrep; ++it) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (*(nim_layers_data + nxy * iz + nx * ix + iy) > 0 && *(nim_columns_data + nxy * iz + nx * ix + iy) > 0) {
                        layer_zi = *(nim_layers_data + nxy * iz + nx * ix + iy)-1;
                        column_yi = *(nim_columns_data + nxy * iz + nx * ix + iy)-1;
                        depth_xi = iz;
                        value_ofinput_data = *(nim_data_data + nxyz * it + nxy * iz + nx * ix + iy) / *(imagiro_vnr_data + nxy_imagiro * layer_zi + nx_imagiro * depth_xi + column_yi);
                        // cout << "value_ofinput_data " << value_ofinput_data << endl;
                        *(imagiro_data + nxyz_imagiro * it + nxy_imagiro * layer_zi + nx_imagiro * depth_xi + column_yi) = *(imagiro_data + nxyz_imagiro * it + nxy_imagiro * layer_zi + nx_imagiro * depth_xi + column_yi) + value_ofinput_data;
                    }
                }
            }
        }
    }
    ///////////////////
    // Fixing wholes //
    ///////////////////

    cout << "  Fixing holes..." << nrep << endl;
    int vinc = 2;  // vicinity
    float value_to_fill = 0;  // Average value in vicinity
    int number_of_vinces = 0;  // Number of non-zero voxels in vicinity

    for (int iz = 0; iz < sizeLayer_imagiro; ++iz) {
        for (int iy = 0; iy < sizeColumn_imagiro; ++iy) {
            for (int ix = 0; ix < sizeDepth_imagiro; ++ix) {
                // cout << iz << " " << ix << " " << iy << endl;
                if (*(imagiro_vnr_data + nxy_imagiro * iz + nx_imagiro * ix + iy) == 0) {
                    // cout << "  iz " << iz << "  ix " << ix << "  iy " << iy << endl;
                    for (int it = 0; it < nrep; ++it) {
                        value_to_fill = 0;
                        number_of_vinces = 0;
                        for (int iy_i=max(0, iy-vinc); iy_i < min(iy+vinc+1, sizeColumn_imagiro); ++iy_i) {
                            for (int ix_i=max(0, ix-vinc); ix_i < min(ix+vinc+1, sizeDepth_imagiro); ++ix_i) {
                                for (int iz_i=max(0, iz-vinc); iz_i < min(iz+vinc+1, sizeLayer_imagiro); ++iz_i) {
                                    if (*(imagiro_vnr_data + nxy_imagiro * iz_i + nx_imagiro * ix_i + iy_i) != 0) {
                                        number_of_vinces = number_of_vinces +1;
                                        value_to_fill = (float) value_to_fill + *(imagiro_data + nxyz_imagiro * it + nxy_imagiro * iz_i + nx_imagiro * ix_i + iy_i);
                                        // cout << "value_to_fill" << value_to_fill << endl;
                                    }
                                }
                            }
                        }
                        *(imagiro_data + nxyz_imagiro * it + nxy_imagiro * iz + nx_imagiro * ix + iy) = value_to_fill/ (float) number_of_vinces;
                    }
                }
            }
        }
    }

    ////////////////////////////////////////
    // Taking care of scale factor if any //
    ////////////////////////////////////////
    cout << "  Writing output 1" << endl;

    // if (nim_data_r->scl_inter == 0){
    imagiro->scl_slope = nim_data_r->scl_slope;
    imagiro->iname_offset = nim_data_r->iname_offset;
    // }

    if (nim_data_r->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    cout << "  Writing output 3 " << endl;

    string prefix = "unfolded_";
    string filename = (string) (data_filename);
    string outfilename = prefix+filename;

    cout << "  Writing as = " << outfilename.c_str() << endl;
    const char * fout_1 = outfilename.c_str();
    if (nifti_set_filenames(imagiro, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(imagiro);

    cout << "  Writing output 2 " << endl;
    const char * fout_5 = "Number_of_voxels.nii";
    if (nifti_set_filenames(imagiro_vnr, fout_5, 1, 1)) {
        return 1;
    }
    nifti_image_write(imagiro_vnr);

    cout << "  Finished." << endl;
    return 0;
}
