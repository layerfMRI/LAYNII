// Note to Renzo. tjere are the comments about copying the same voxel to itself when the distance is smaller than 1.7. 
// remove this stuff

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_COLUMNAR_DIST : Calculates cortical distances (columnar structures) \n"
    "                   based on the gray matter geometry.\n"
    "\n"
    "Usage:\n"
    "    LN_COLUMNAR_DIST -layers layers.nii -landmarks landmarks.nii \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -layers     : Nifti (.nii) that contains layer or column masks. \n"
    "    -landmarks  : Nifti (.nii) with landmarks (use value 1 as origin). \n"
    "                  Landmarks should be at least 4 voxels thick. \n"
    "    -vinc       : Maximal length of cortical distance. Bigger values \n"
    "                  will take longer but span a larger cortical area. \n"
    "                  Default is 40. \n"
    "    -Ncolumns   : (Optional) For the number of columns. Smaller values \n"
    "                  will result in thick columns \n"
    "    -verbose    : (Optional) to write out all the intermediate \n"
    "                  steps of the algorithm (e.g. for debugging) \n"
    "\n"
    "Notes:\n"
    "    - The layer nii file and the landmarks nii file should have the \n"
    "      same dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
    char *fin_layer = NULL, *fin_landmark = NULL;
    int ac, vinc_max = 40, Ncolumns = 0;
    int verbose = 0;
    if (argc < 3) return show_help();

    // process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -layers\n");
                return 1;
            }
            fin_layer = argv[ac];
        } else if (!strcmp(argv[ac], "-landmarks")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -landmarks\n");
                return 1;
            }
            fin_landmark = argv[ac];
        } else if (!strcmp(argv[ac], "-vinc")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -vinc\n");
                return 1;
            }
            vinc_max = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-Ncolumns")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Ncolumns\n");
                return 1;
            }
            Ncolumns = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-verbose")) {
            verbose = 1;
            cout << "  Debug mode active. More outputs." << endl;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_layer) {
        fprintf(stderr, "** missing option '-layers'\n");
        return 1;
    }
    if (!fin_landmark) {
        fprintf(stderr, "** missing option '-landmarks'\n");
        return 1;
    }

    // Read input dataset
    nifti_image * nii_input1 = nifti_image_read(fin_layer, 1);
    if (!nii_input1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_layer);
        return 2;
    }
    nifti_image * nii_input2 = nifti_image_read(fin_landmark, 1);
    if (!nii_input2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_landmark);
        return 2;
    }

    log_welcome("LN_COLUMNAR_DIST");
    log_nifti_descriptives(nii_input1);  // Layers
    log_nifti_descriptives(nii_input2);  // Landmarks

    // Get dimensions of input
    const int size_x = nii_input1->nx;
    const int size_y = nii_input1->ny;
    const int size_z = nii_input1->nz;
    const int nr_voxels = size_z * size_y * size_x;
    const int nx = nii_input1->nx;
    const int nxy = nii_input1->nx * nii_input1->ny;
    const int nxyz = nii_input1->nx * nii_input1->ny * nii_input1->nz;
    const float dX = nii_input1->pixdim[1];
    const float dY = nii_input1->pixdim[2];
    const float dZ = nii_input1->pixdim[3];

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii_layers = copy_nifti_as_float32(nii_input1);
    float* nim_layers_data = static_cast<float*>(nii_layers->data);
    nifti_image* nii_landmarks = copy_nifti_as_int32(nii_input2);
    int* nim_landmarks_data = static_cast<int*>(nii_landmarks->data);

    // Finding number of layers
    int nr_layers = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_layers_data + i) > nr_layers) {
            nr_layers = *(nim_layers_data + i);
        }
    }
    cout << "  There are " << nr_layers << " layers." << endl;

    // Allocating necessary files
    nifti_image* Grow_x = copy_nifti_as_int32(nii_layers);
    nifti_image* Grow_y = copy_nifti_as_int32(nii_layers);
    nifti_image* Grow_z = copy_nifti_as_int32(nii_layers);

    int32_t* Grow_x_data = static_cast<int32_t*>(Grow_x->data);
    int32_t* Grow_y_data = static_cast<int32_t*>(Grow_y->data);
    int32_t* Grow_z_data = static_cast<int32_t*>(Grow_z->data);

    nifti_image * growfromCenter = copy_nifti_as_int32(nii_layers);
    nifti_image * growfromCenter_thick = copy_nifti_as_int32(nii_layers);

    int32_t* growfromCenter_data = static_cast<int32_t*>(growfromCenter->data);
    int32_t* growfromCenter_thick_data = static_cast<int32_t*>(growfromCenter_thick->data);

    for (int i = 0; i < nr_voxels; ++i) {
        *(Grow_x_data + i) = 0;
        *(Grow_y_data + i) = 0;
        *(Grow_z_data + i) = 0;
        *(growfromCenter_data + i) = 0;
        *(growfromCenter_thick_data + i) = 0;
    }

    // ========================================================================
    // Prepare growing variables
    // ========================================================================
    float x1g = 0., y1g = 0., z1g = 0.;

    float min_val = 0.;
    float dist_min2 = 0.;
    float dist_i = 0.;
    float dist_p1;

    int grow_vinc_area = 1;

    // ========================================================================
    // Growing from Center cross columns
    // ========================================================================
    int jz_start, jy_start, jx_start, jz_stop, jy_stop, jx_stop;
    int kz_start, ky_start, kx_start, kz_stop, ky_stop, kx_stop;

    cout << "  Growing from center..." << endl;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_y; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                // Defining seed at center landmark
                if (*(nim_landmarks_data + voxel_i) == 1
                    && abs((int) (*(nim_layers_data + voxel_i) - nr_layers / 2)) < 2) {
                    *(growfromCenter_data + voxel_i) = 1.;
                    *(Grow_x_data + voxel_i) = ix;
                    *(Grow_y_data + voxel_i) = iy;
                    *(Grow_z_data + voxel_i) = iz;
                }
            }
        }
    }

    for (int grow_i = 1; grow_i < vinc_max; grow_i++) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxy * iz + nx * iy + ix;

                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;

                    if (abs((int) (*(nim_layers_data + voxel_i) - nr_layers / 2)) < 2
                        && *(growfromCenter_data + voxel_i) == 0
                        && *(nim_landmarks_data + voxel_i) < 2) {
                        // Only grow into areas that are GM and that have not
                        // been grown into, yet... and it should stop as soon
                        // as it hits the border

                        jz_start = max(0, iz - grow_vinc_area);
                        jz_stop = min(iz + grow_vinc_area, size_z - 1);
                        jy_start = max(0, iy - grow_vinc_area);
                        jy_stop = min(iy + grow_vinc_area, size_y - 1);
                        jx_start = max(0, ix - grow_vinc_area);
                        jx_stop = min(ix + grow_vinc_area, size_x - 1);

                        for (int jz = jz_start; jz <= jz_stop; ++jz) {
                            for (int jy = jy_start; jy <= jy_stop; ++jy) {
                                for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                    int voxel_j = nxy * jz + nx * jy + jx;

                                    dist_i = dist((float)ix, (float)iy, (float)iz,
                                                  (float)jx, (float)jy, (float)jz,
                                                  dX, dY, dZ);

                                    if (*(growfromCenter_data + voxel_j) == grow_i
                                        && *(nim_landmarks_data + voxel_i) < 2) {
                                        if (dist_i < dist_min2) {
                                            dist_min2 = dist_i;
                                            x1g = jx;
                                            y1g = jy;
                                            z1g = jz;
                                            dist_p1 = dist_min2;
                                        }
                                    }
                                }
                            }
                        }
                        if (dist_min2 < 1.7) {  // ???? I DONT REMEMBER WHY I NEED THIS ????
                            *(growfromCenter_data + voxel_i) = grow_i+1;
                            *(Grow_x_data + voxel_i) = *(Grow_x_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_y_data + voxel_i) = *(Grow_y_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                            *(Grow_z_data + voxel_i) = *(Grow_z_data + nxy * (int)z1g + nx * (int)x1g + (int)y1g);
                        }
                    }
                }
            }
        }
    }
    if (verbose == 1) {
        save_output_nifti(fin_layer, "coordinates_1_path", growfromCenter, false);
    }

    // ========================================================================
    // Smooth columns
    // ========================================================================
    // NOTE(Renzo): In the future this should be done only within connected
    // areas. Right now there might be a problem, when the center of two GM.
    cout << "  Smoothing in middle layer..." << endl;

    nifti_image* smooth = copy_nifti_as_float32(nii_layers);
    float* smooth_data = static_cast<float*>(smooth->data);

    int FWHM_val = 1;
    int vinc_sm = 5;  // if voxel is too far away, I ignore it.
    dist_i = 0.;
    cout << "    vinc_sm " << vinc_sm << endl;
    cout << "    FWHM_val " << FWHM_val << endl;
    cout << "    Here 2 " << endl;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                float total_weight = 0;
                *(smooth_data + voxel_i) = 0 ;

                if (*(growfromCenter_data + voxel_i)  > 0) {

                    jz_start = max(0, iz - vinc_sm);
                    jy_start = max(0, iy - vinc_sm);
                    jx_start = max(0, ix - vinc_sm);
                    jz_stop = min(iz + vinc_sm, size_z - 1);
                    jy_stop = min(iy + vinc_sm, size_y - 1);
                    jx_stop = min(ix + vinc_sm, size_x - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;

                                if (*(growfromCenter_data + voxel_j)  > 0) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz,
                                                  (float)jx, (float)jy, (float)jz,
                                                  dX, dY, dZ);

                                    float w = gaus(dist_i, FWHM_val);
                                    *(smooth_data + voxel_i) += *(growfromCenter_data + voxel_j) * w;
                                    total_weight += w;
                                }
                            }
                        }
                    }
                    if (total_weight > 0) {
                        *(smooth_data + voxel_i) /= total_weight;
                    }
                }
            }
        }
    }

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(growfromCenter_data + i) > 0) {
            *(growfromCenter_data + i) = static_cast<int>(*(smooth_data + i));
        }
    }

    if (verbose == 1) {
        save_output_nifti(fin_layer, "coordinates_2_path_smooth",
                          growfromCenter, false);
    }

    // ========================================================================
    // Extending columns across layers
    // ========================================================================
    // NOTE(Renzo): This is not perfect yet, because it has only 4 directions
    // to grow thus ther might be orientation biases.
    cout << "  Extending columns across layers..." << endl;

    nifti_image* hairy = copy_nifti_as_int32(nii_layers);
    int32_t* hairy_data = static_cast<int32_t*>(hairy->data);

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_layers_data + i) > 1) {
            *(hairy_data + i) = 1;
        } else {
            *(hairy_data + i) = 0;
        }
    }

    // This is an upper limit of the cortical thickness
    dist_min2 = 10000.;

    // This is step size that neigbouring GM voxels need to be to be classified
    // as one side of the GM bank.
    int vinc_steps = 1;

    int vinc_sm_g = 25;
    int pref_ratio = 0;

    // For estimation of time
    int nvoxels_to_go_across = 0;
    int running_index = 0;
    for (int i = 0; i < nr_voxels; ++i) {
        if (*(nim_layers_data + i) > 1
            && *(nim_layers_data + i) < nr_layers - 1) {
                nvoxels_to_go_across++;
        }
    }

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                if (*(nim_layers_data + voxel_i) > 1
                    && *(nim_layers_data + voxel_i) < nr_layers - 1) {
                    running_index++;
                    int temp = (running_index * 100) / nvoxels_to_go_across;
                    if (temp != pref_ratio) {
                        cout << "\r    " << temp << "% " << flush;
                        pref_ratio = temp;
                    }
                    // --------------------------------------------------------
                    // Find area that is not from the other sulcus
                    // --------------------------------------------------------
                    jz_start = max(0, iz - vinc_sm_g - vinc_steps);
                    jy_start = max(0, iy - vinc_sm_g - vinc_steps);
                    jx_start = max(0, ix - vinc_sm_g - vinc_steps);
                    jz_stop = min(iz + vinc_sm_g + vinc_steps, size_z - 1);
                    jy_stop = min(iy + vinc_sm_g + vinc_steps, size_y - 1);
                    jx_stop = min(ix + vinc_sm_g + vinc_steps, size_x - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;
                                *(hairy_data + voxel_j) = 0;
                            }
                        }
                    }
                    min_val = 0;

                    // Iteration loop that determines a local patch of connected
                    // voxels, excluding voxels from opposite GM bank.
                    // NOTE(Renzo): This loop takes forever
                    *(hairy_data + voxel_i) = 1;
                    for (int K_ = 0; K_ < vinc_sm_g; K_++) {

                        jz_start = max(0, iz - vinc_sm_g);
                        jy_start = max(0, iy - vinc_sm_g);
                        jx_start = max(0, ix - vinc_sm_g);
                        jz_stop = min(iz + vinc_sm_g, size_z - 1);
                        jy_stop = min(iy + vinc_sm_g, size_x - 1);
                        jx_stop = min(ix + vinc_sm_g, size_y - 1);

                        for (int jz = jz_start; jz <= jz_stop; ++jz) {
                            for (int jy = jy_start; jy <= jy_stop; ++jy) {
                                for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                    int voxel_j = nxy * jz + nx * jy + jx;

                                    if (*(hairy_data + voxel_j) == 1) {

                                        kz_start = max(0, jz - vinc_steps);
                                        ky_start = max(0, jy - vinc_steps);
                                        kx_start = max(0, jx - vinc_steps);
                                        kz_stop = min(jz + vinc_steps, size_z - 1);
                                        ky_stop = min(jy + vinc_steps, size_y - 1);
                                        kx_stop = min(jx + vinc_steps, size_x - 1);

                                        for (int kz = kz_start; kz <= kz_stop; ++kz) {
                                            for (int ky = ky_start; ky <= ky_stop; ++ky) {
                                                for (int kx = kx_start; kx <= kx_stop; ++kx) {
                                                    int voxel_k = nxy * kz + nx * ky + kx;

                                                    float d = dist((float)jx, (float)jy, (float)jz,
                                                                   (float)kx, (float)ky, (float)kz,
                                                                   1, 1, 1);

                                                    if (d <= 1
                                                        && *(nim_layers_data + voxel_k) > 1
                                                        && *(nim_layers_data + voxel_k) < nr_layers - 1) {
                                                        *(hairy_data + voxel_k) = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    // Only grow into areas that are GM and that have not been
                    // grown into, yet... and it should stop as soon as it hits
                    // the border

                    jz_start = max(0, iz - vinc_sm_g);
                    jy_start = max(0, iy - vinc_sm_g);
                    jx_start = max(0, ix - vinc_sm_g);
                    jz_stop = min(iz + vinc_sm_g, size_z - 1);
                    jy_stop = min(iy + vinc_sm_g, size_x - 1);
                    jx_stop = min(ix + vinc_sm_g, size_y - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;

                                dist_i = dist((float)ix, (float)iy, (float)iz,
                                              (float)jx, (float)jy, (float)jz,
                                              dX, dY, dZ);

                                if (*(hairy_data + voxel_j) == 1
                                    && dist_i < dist_min2
                                    && *(growfromCenter_data + voxel_j) > 0) {
                                    dist_min2 = dist_i;
                                    x1g = jx;
                                    y1g = jy;
                                    z1g = jz;
                                    dist_p1 = dist_min2;
                                    min_val = *(growfromCenter_data + voxel_j);
                                }
                            }
                        }
                    }
                    *(growfromCenter_thick_data + voxel_i) = min_val;
                }
            }
        }
    }
    cout << endl;  // to close the online output

    if (verbose == 1) {
        save_output_nifti(fin_layer, "coordinates_3_thick",
                          growfromCenter_thick, false);
    }

    ////////////////////////////////////////////////////////////
    // Smooth columns within the thick cortex.                //
    // This smoothing is done to correct for Pytagoras errors //
    // This smoothing is within GM banks only                 //
    // (which makes it slow)                                  //
    ////////////////////////////////////////////////////////////
    cout << "  Smoothing the thick cortex now..." << endl;
    for (int i = 0; i < nr_voxels; ++i) {
        *(smooth_data + i) = 0;
    }

    FWHM_val = 3;

    // This is small so that neighboring sulci are not affecting each other
    vinc_sm = 8;  // max(1., 2. * FWHM_val/dX);

    cout << "    vinc_sm " << vinc_sm << endl;
    cout << "    FWHM_val " << FWHM_val << endl;
    cout << "    Starting extended now  " << endl;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                float total_weight = 0;
                // *(smooth_data + voxel_i) = 0;

                if (*(growfromCenter_thick_data + voxel_i) > 0) {
                    // --------------------------------------------------------
                    // Find area that is not from the other sulcus
                    // --------------------------------------------------------
                    // Preparation of dummy vicinity
                    jz_start = max(0, iz - vinc_sm - vinc_steps);
                    jy_start = max(0, iy - vinc_sm - vinc_steps);
                    jx_start = max(0, ix - vinc_sm - vinc_steps);
                    jz_stop = min(iz + vinc_sm + vinc_steps, size_z - 1);
                    jy_stop = min(iy + vinc_sm + vinc_steps, size_y - 1);
                    jx_stop = min(ix + vinc_sm + vinc_steps, size_x - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;
                                *(hairy_data + voxel_j) = 0;
                            }
                        }
                    }
                    *(hairy_data + voxel_i) = 1;
                    for (int K_ = 0; K_ < vinc_sm; K_++) {

                        jz_start = max(0, iz - vinc_sm);
                        jy_start = max(0, iy - vinc_sm);
                        jx_start = max(0, ix - vinc_sm);
                        jz_stop = min(iz + vinc_sm, size_z - 1);
                        jy_stop = min(iy + vinc_sm, size_y - 1);
                        jx_stop = min(ix + vinc_sm, size_x - 1);

                        for (int jz = jz_start; jz <= jz_stop; ++jz) {
                            for (int jy = jy_start; jy <= jy_stop; ++jy) {
                                for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                    int voxel_j = nxy * jz + nx * jy + jx;

                                    if (*(hairy_data + voxel_j) == 1) {

                                        int kz_start = max(0, jz - vinc_steps);
                                        int ky_start = max(0, jy - vinc_steps);
                                        int kx_start = max(0, jx - vinc_steps);
                                        int kz_stop = min(jz + vinc_steps, size_z - 1);
                                        int ky_stop = min(jy + vinc_steps, size_y - 1);
                                        int kx_stop = min(jx + vinc_steps, size_x - 1);

                                        for (int kz = kz_start; kz <= kz_stop; ++kz) {
                                            for (int ky = ky_start; ky <= ky_stop; ++ky) {
                                                for (int kx = kx_start; kx <= kx_stop; ++kx) {
                                                    int voxel_k = nxy * kz + nx * ky + kx;

                                                    float d = dist((float)jx, (float)jy, (float)jz,
                                                                   (float)kx, (float)ky, (float)kz,
                                                                   1, 1, 1);
                                                    if (d <= 1
                                                        && *(nim_layers_data + voxel_k) > 1
                                                        && *(nim_layers_data + voxel_k) < nr_layers - 1) {
                                                        *(hairy_data + voxel_k) = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Smoothing within each layer and within the local patch
                    int nr_layers_i = *(nim_layers_data + voxel_i);

                    jz_start = max(0, iz - vinc_sm);
                    jy_start = max(0, iy - vinc_sm);
                    jx_start = max(0, ix - vinc_sm);
                    jz_stop = min(iz + vinc_sm, size_z - 1);
                    jy_stop = min(iy + vinc_sm, size_y - 1);
                    jx_stop = min(ix + vinc_sm, size_x - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;

                                if (*(hairy_data + voxel_j) == 1
                                    && abs((int) *(nim_layers_data + voxel_j) - nr_layers_i) < 2
                                    && *(growfromCenter_thick_data + voxel_j) > 0) {
                                    dist_i = dist((float)ix, (float)iy, (float)iz,
                                                  (float)jx, (float)jy, (float)jz,
                                                  dX, dY, dZ);
                                    float w = gaus(dist_i, FWHM_val);
                                    *(smooth_data + voxel_i) += *(growfromCenter_thick_data + voxel_j) * w;
                                    total_weight += w;
                                }
                            }
                        }
                    }
                    if (total_weight > 0) {
                        *(smooth_data + voxel_i) /= total_weight;
                    }
                }
            }
        }
    }
    cout << "    Extended now." << endl;

    for (int i = 0; i < nr_voxels; ++i) {
        if (*(growfromCenter_thick_data + i) > 0) {
            *(growfromCenter_thick_data + i) = *(smooth_data + i);
        }
    }
    cout << "    Smoothing done." << endl;

    if (verbose == 1) {
        save_output_nifti(fin_layer, "coordinates_4_thick_smooth",
                          growfromCenter_thick, false);
    }

    // ========================================================================
    // Grow final outer rim
    // ========================================================================
    // Note(Renzo): I could not fill this earlier because it would have resulted
    // in leakage from the opposite back. Thus I always left a safety corridor.
    // Which is filled in now.
    int vinc_rim = 2;
    for (int i = 0; i < nr_voxels; ++i) {
        *(hairy_data + i) = *(growfromCenter_thick_data + i);
    }

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;

                if (*(nim_layers_data + voxel_i) > 0
                    && *(growfromCenter_thick_data + voxel_i) == 0) {
                    dist_min2 = 10000.;
                    x1g = 0;
                    y1g = 0;
                    z1g = 0;
                    min_val = 0;
                    // Only grow into areas that are GM and that have not been
                    // grown into, yet... And it should stop as soon as it hits
                    // the border.

                    jz_start = max(0, iz - vinc_rim);
                    jy_start = max(0, iy - vinc_rim);
                    jx_start = max(0, ix - vinc_rim);
                    jz_stop = min(iz + vinc_rim, size_z - 1);
                    jy_stop = min(iy + vinc_rim, size_y - 1);
                    jx_stop = min(ix + vinc_rim, size_x - 1);

                    for (int jz = jz_start; jz <= jz_stop; ++jz) {
                        for (int jy = jy_start; jy <= jy_stop; ++jy) {
                            for (int jx = jx_start; jx <= jx_stop; ++jx) {
                                int voxel_j = nxy * jz + nx * jy + jx;

                                dist_i = dist((float)ix, (float)iy, (float)iz,
                                              (float)jx, (float)jy, (float)jz,
                                              dX, dY, dZ);

                                if (dist_i < dist_min2
                                    && *(nim_layers_data + voxel_j) > 0
                                    && *(growfromCenter_thick_data + voxel_j) > 0) {
                                    dist_min2 = dist_i;
                                    x1g = jx;
                                    y1g = jy;
                                    z1g = jz;
                                    dist_p1 = dist_min2;
                                    min_val = *(growfromCenter_thick_data + voxel_j);
                                }
                            }
                        }
                    }
                    *(hairy_data + voxel_i) = min_val;
                }
            }
        }
    }

    if (verbose == 1) {
        save_output_nifti(fin_layer, "coordinates_5_extended", hairy, false);
    }

    // ========================================================================
    // Resample the number of columns to the user given value
    // ========================================================================
    if (Ncolumns > 0)  {
        cout << "   Resampling the number of columns..." << endl;
        int max_columns = 0;
        int min_columns = 100000000;
        for (int i = 0; i < nr_voxels; ++i) {
            int val = static_cast<int>(*(hairy_data + i));
            if (val >  0) {
                if (val >  max_columns) {
                    max_columns = val;
                }
                if (val < min_columns) {
                    min_columns = val;
                }
            }
        }

        cout << "   Max = " << max_columns << " | Min = " << min_columns << endl;
        for (int i = 0; i < nr_voxels; ++i) {
            if (*(hairy_data + i) >  0) {
                int a = static_cast<int>(*(hairy_data + i) - min_columns);
                int b = static_cast<int>(Ncolumns - 1);
                int c = static_cast<int>(max_columns - min_columns);
                *(hairy_data + i) = a * b / c + 1;
            }
        }
    }
    save_output_nifti(fin_layer, "coordinates_final", hairy, true);

    cout << "  Finished." << endl;
    return 0;
}
