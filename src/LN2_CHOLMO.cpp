
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_CHOLMO: Add additional layers onto already existing layers with\n"
    "            like the name suggests is is padding layers \n"
    "            this can be helpful to extend profiles into CSF or into WM \n"
    "\n"
    "Usage:\n"
    "    LN2_CHOLMO -layers layers.nii \n"
    "    LN2_CHOLMO -layers layers.nii -outer -nr_layers 3 -layer_thickness 0.5 \n"
    "    LN2_CHOLMO -layers layers.nii -inner -nr_layers 2 -layer_thickness 0.4\n"
    "    ../LN2_CHOLMO -layers sc_layers.nii -outer -nr_layers 3 -layer_thickness 0.4 -output padded_layers.nii \n"
    "\n"
    "Options:\n"
    "    -help            : Show this help.\n"
    "    -layers          : Specify input dataset of layers\n"
    "                       It is assumed that it conists of intager numbers of layers\n"
    "                       It is assumed that deeper layers have small values \n"
    "                       It is assumes that superficial layers have large values.\n"
    "    -nr_layers       : Number of layers to add on top of the already existing ones.\n"
    "                       Default is 3.\n"
    "    -layer_thickness : Thickness of the added layers in mm. Default is 0.8.\n"
    "                       This should not be smaller than the voxel dimension.\n"
    "    -debug           : (Optional) Save extra intermediate outputs.\n"
    "    -output          : (Optional) Output basename. Default is '_padded' as suffix.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin = NULL, *fout = NULL;
    uint16_t ac, nr_layers = 3;
    double_t layer_thickness = 0.8;
    bool  mode_debug = false,  mode_inner = false,  mode_outer = true;
    bool  use_outpath = false;


    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layers\n");
                return 1;
            }
            fin = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-nr_layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_layers\n");
            } else {
                nr_layers = atof(argv[ac]);
            }
       } else if (!strcmp(argv[ac], "-layer_thickness")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -interations\n");
            } else {
                layer_thickness = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-inner")) {
            mode_inner = true;
            mode_outer = false;
        } else if (!strcmp(argv[ac], "-outer")) {
            mode_outer = true;
            mode_inner = false;
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-layers'\n");
        return 1;
    }

    if (mode_outer == mode_inner) {
        fprintf(stderr, "** You selected to pad the layers both at the inner \n"
                        "   and at the outer GM borders?!?!!?? \n"
                        "   Make up you mind!\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN2_CHOLMO");
    log_nifti_descriptives(nii1);

    cout << "  Nr. layers that should be added: " << nr_layers << endl;
    cout << "  Thickness of those added layers: "  << layer_thickness << endl;


    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_layers = copy_nifti_as_int16(nii1);
    int16_t* nii_layers_data = static_cast<int16_t*>(nii_layers->data);

    nifti_image* step = copy_nifti_as_int16(nii_layers);
    int16_t* step_data = static_cast<int16_t*>(step->data);

    nifti_image* dist = copy_nifti_as_float32(nii_layers);
    float* dist_data = static_cast<float*>(dist->data);


    // ========================================================================
    // Looking what is already there in input
    // ========================================================================
    int max_layers = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_layers_data + i) >= max_layers){
            max_layers = *(nii_layers_data + i);
        }
    }
    cout << "  There are " << max_layers<< " layers already. I will add another  " << nr_layers << endl;

    // ========================================================================
    // Grow outwards
    // ========================================================================
    cout << "\n  Start growing ....." << endl;

    if (mode_outer) {
        // Initialize grow volume for outwards growign
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_layers_data + i) == max_layers) {  // starting point to grow at the outer most layer with the largest number
                *(step_data + i) = 1;
                *(dist_data + i) = 0.;
            } else {
                *(step_data + i) = 0.;
                *(dist_data + i) = 0.;
            }
        }
    }


    if (mode_inner) {
        // Initialize grow volume for outwards growign
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_layers_data + i) == 1) {  // starting point to grow at the inner most layer with value 1
                *(step_data + i) = 1;
                *(dist_data + i) = 0.;
            } else {
                *(step_data + i) = 0.;
                *(dist_data + i) = 0.;
            }
        }
    }

    uint16_t grow_step = 1;
    uint32_t voxel_counter = nr_voxels;
    uint32_t ix, iy, iz, j, k;
    float d;
    while (voxel_counter != 0) {
        cout << "\r  Growing step " << grow_step << "......"  << flush;
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {

            if (*(step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours (faces of cube)
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dX;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dX;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dY;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dY;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dZ;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dZ;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours (edges of cube)
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xy;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xy;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xy;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xy;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_yz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_yz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_yz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_yz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours (corners of cube)
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j)  || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_layers_data + j) == 0) {
                        d = *(dist_data + i) + dia_xyz;
                        if (d < *(dist_data + j) || *(dist_data + j) == 0) {
                            *(dist_data + j) = d;
                            *(step_data + j) = grow_step + 1;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }

    if (mode_debug) {
        save_output_nifti(fout, "step", step, false);
        save_output_nifti(fout, "dist", dist, false);
    }

    // ========================================================================
    // Quantizing the grown distances into layers of desired thickness
    // ========================================================================

    for (uint32_t i = 0; i != nr_voxels; ++i) {
            // truncate the metric to layers
            *(dist_data + i) = ceil( *(dist_data + i) / layer_thickness) ;
    }

    // ========================================================================
    // Giving some feedback to the user that might help with interpretation.
    // ========================================================================
    int max_layers_padded = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(dist_data + i) >= max_layers_padded){
            max_layers_padded = *(dist_data + i);
        }
    }
    if (max_layers_padded < nr_layers ){
        cout << "\n  Comments about the output:" << endl;
        cout << "    -> there was not enough space to add as many layers " << endl;
        cout << "       with the desired thickness.  " << endl;
        cout << "       There was only space for " << max_layers_padded<< " additional layers, but you wanted "<<nr_layers << " additional layers" << endl;
        cout << "       This can be due to the fact that there brain fills the entire FOV for outwards growing.  " << endl;
        cout << "       This can be due to the fact that there brain fills the WM for inwards growing.  " << endl;
        cout << "       Consider thinner layers, or less layers  " << endl;
        cout << "       Consider adding empty slces to increase the FOV: ImageMath 3 padded.nii PadImage data.nii 200 0" << endl;
    }

   if (layer_thickness < dX || layer_thickness < dY || layer_thickness < dZ ){
        cout << "\n  Comments about the output:" << endl;
        cout << "    -> The layers might be too thin for the voxel grid of the data." << endl;
        cout << "       You chose a layer thickness of " << layer_thickness << " mm " <<endl;
        cout << "       while the voxel grid is (X,Y,Z) = (" <<  dX << "," << dY << "," << dZ<< ")" << endl;
        cout << "       Consider using thicker layers  " << endl;
        cout << "       Consider upsampling the layer file: " << endl;
        cout << "           https://github.com/ofgulban/LAYNII_extras/blob/master/Upsamp_2rimify/Up_sample_2d.sh " << endl;
   }

    // ========================================================================
    // Truncate the layers that are not desired
    // ========================================================================
    //Writing out un-trunkated layers
    if (mode_debug) {
        save_output_nifti(fout, "dist_quantized", dist, false);
    }

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(dist_data + i) >  nr_layers) {
            *(dist_data + i) = 0.0 ;
        }
    }

    // ========================================================================
    // Add desired layers to what is already there for outwards growing
    // ========================================================================
    cout << "\n  combining already existing layers with the padded ones..." << endl << endl;

    if (mode_outer) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(dist_data + i) != 0 ) {
                *(dist_data + i) =  *(dist_data + i)  + max_layers;
                } else {
                *(dist_data + i) = *(nii_layers_data + i) ;
                }
        }
    }

    if (mode_inner) {
        if ( nr_layers > max_layers_padded ) {
            cout << "\n  Since there is not enough space to add " << nr_layers << "layers (see above)," << endl;
            cout << "  I am just using " << max_layers_padded << "layers instead" << endl << endl;
                    nr_layers = max_layers_padded ;
        }

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(dist_data + i) != 0 ) {
                *(dist_data + i) = nr_layers - *(dist_data + i) + 1 ;
            }
            if (*(dist_data + i) == 0 && *(nii_layers_data + i) != 0 ) {
                *(dist_data + i) = *(nii_layers_data + i) + nr_layers  ;
            }
        }
    }

    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "padded", dist, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
