#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_MULTILATERATE: Injects a coordinate system upon a region of the rim file.\n"
    "                   These coordinates can be used to flatten chunks of the brain.\n"
    "                   Or, to generate bins/cells (subsets of voxels).\n"
    "\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "!!!!    WORK IN PROGRESS    !!!!\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "\n"
    "Usage:\n"
    "    LN2_MULTILATERATE -rim rim.nii -midgm rim_midgm.nii -centroid rim_midgm_centroid.nii -radius 10\n"
    "    LN2_MULTILATERATE -rim rim.nii -midgm rim_midgm.nii -centroid rim_midgm_extrema.nii\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -rim      : Segmentation input. Use 3 to code gray matter voxels\n"
    "                This program only injects coordinates to the voxels\n"
    "                labeled with 3.\n"
    "    -midgm    : Middle gray matter file (from LN2_LAYERS output).\n"
    "    -centroid : A nifti file that contains either of the following:\n"
    "                (I) One voxel labeled with value '1' for indicating the\n"
    "                centroid of the region of interest; or (II) four voxels\n"
    "                each labeled with values 1, 2, 3, 4. In case (I), you\n"
    "                need to use '-radius' parameter. In case (II), voxels\n"
    "                labeled with 1 & 2 determine the first axis of the\n"
    "                coordinate system and the voxels labeled with 3 & 4\n"
    "                determine the second axis."
    "    -radius   : Distance threshold from centroid.\n"
    "    -debug    : (Optional) Save extra intermediate outputs.\n"
    "    -output   : (Optional) Output basename for all outputs.\n"
    "\n"
    "Notes:\n"
    "    - This program is written for 3D images. We might add 2D image support\n"
    "      in the future."
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL;
    float thr_radius = 10;
    int ac;
    bool mode_debug = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -rim\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-midgm")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -midgm\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-centroid")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -centroid\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-radius")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radius\n");
                return 1;
            }
            thr_radius = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-midgm'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-centroid'\n");
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
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }

    log_welcome("LN2_MULTILATERATE");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

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
    nifti_image* nii_rim = copy_nifti_as_int32(nii1);
    int32_t* nii_rim_data = static_cast<int32_t*>(nii_rim->data);
    nifti_image* nii_midgm = copy_nifti_as_int32(nii2);
    int32_t* nii_midgm_data = static_cast<int32_t*>(nii_midgm->data);
    nifti_image* nii_centroid = copy_nifti_as_int32(nii3);
    int32_t* nii_centroid_data = static_cast<int32_t*>(nii_centroid->data);

    // Prepare required nifti images
    nifti_image* quartet = copy_nifti_as_int32(nii_rim);
    int32_t* quartet_data = static_cast<int32_t*>(quartet->data);
    nifti_image* flood_step = copy_nifti_as_int32(nii_rim);
    int32_t* flood_step_data = static_cast<int32_t*>(flood_step->data);
    nifti_image* flood_dist = copy_nifti_as_float32(nii_rim);
    float* flood_dist_data = static_cast<float*>(flood_dist->data);

    nifti_image* perimeter = copy_nifti_as_int32(nii_rim);
    int32_t* perimeter_data = static_cast<int32_t*>(perimeter->data);

    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(quartet_data + i) = 0;
        *(flood_step_data + i) = 0;
        *(flood_dist_data + i) = 0;
        *(perimeter_data + i) = 0;
    }

    // Final voronoi volume to output midgm distances for whole rim
    nifti_image* voronoi = copy_nifti_as_float32(nii_rim);
    float* voronoi_data = static_cast<float*>(voronoi->data);
    nifti_image* smooth = copy_nifti_as_float32(voronoi);
    float* smooth_data = static_cast<float*>(smooth->data);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // --------------------------------------------------------------------
    // NOTE(Faruk) This part is for speeding up the final Voronoi propagation
    // --------------------------------------------------------------------
    // Reduce number of looped-through voxels for the second stage (Voronoi)
    uint32_t nr_voi2 = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3){
            nr_voi2 += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id2;
    voi_id2 = (int32_t*) malloc(nr_voi2*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t iii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3){
            *(voi_id2 + iii) = i;
            iii += 1;
        }
    }

    // ========================================================================
    // Initial flood from centroid
    // ========================================================================
    cout << "  Start trilaterating..." << endl;
    // Find the initial voxel
    uint32_t start_voxel;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_centroid_data + i) == 1) {
            start_voxel = i;
        }
    }
    *(nii_midgm_data + start_voxel) = 2;

    // Find distances from input centroid
    float flood_dist_thr = std::numeric_limits<float>::infinity();
    int32_t grow_step = 1;
    uint32_t voxel_counter = nr_voxels;
    uint32_t ix, iy, iz, i, j;
    float d;

    // Initialize grow volume
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 2) {
            *(flood_step_data + i) = 1.;
            *(flood_dist_data + i) = 0.;
        } else if (*(flood_dist_data + i) >= flood_dist_thr
                   && *(flood_dist_data + i) > 0) {
            *(flood_step_data + i) = 0.;
            *(flood_dist_data + i) = 0.;
            *(nii_midgm_data + i) = 1;
        } else if (*(flood_dist_data + i) < flood_dist_thr
                   && *(flood_dist_data + i) > 0) {
            *(nii_midgm_data + i) = 0;  // no need to recompute
        }
    }

    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);  // Map subset to full set
            if (*(flood_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // --------------------------------------------------------
                // 1-jump neighbours
                // --------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                // --------------------------------------------------------
                // 2-jump neighbours
                // --------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }

                // --------------------------------------------------------
                // 3-jump neighbours
                // --------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_midgm_data + j) == 1) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }

    save_output_nifti(fout, "centroid_dist", flood_dist, true);

    // Translate 0 crossing
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);
        *(flood_dist_data + i) -= thr_radius;
    }

    // ========================================================================
    // Find perimeter
    // ========================================================================
    cout << "\n  Start finding perimeter..." << endl;

    // Mark voxels inside perimeter
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);
        if (*(flood_dist_data + i) < 0) {
            *(perimeter_data + i) = 1;
        }
    }

    // TODO(Faruk): Once the proof of concept is working, I need to reduce code
    // repetition here through functions. Signbit checks are easy reductions.
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);  // Map subset to full set
        // Check sign changes to find zero crossings
        if (*(flood_dist_data + i) == 0) {
            *(perimeter_data + i) = 2;
        } else {
            float m = *(flood_dist_data + i);
            float n;
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

            // --------------------------------------------------------
            // 1-jump neighbours
            // --------------------------------------------------------
            if (ix > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            // --------------------------------------------------------
            // 2-jump neighbours
            // --------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            // --------------------------------------------------------
            // 3-jump neighbours
            // --------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                if (*(nii_midgm_data + j) == 1){
                    n = *(flood_dist_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(perimeter_data + i) = 2;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(perimeter_data + j) = 2;
                        } else {  // Equal +/- normalized distance
                            *(perimeter_data + i) = 2;
                            *(perimeter_data + j) = 2;
                        }
                    }
                }
            }
        }
    }
    save_output_nifti(fout, "perimeter", perimeter, false);

    // --------------------------------------------------------------------
    // Grow perimeter within rim to generate vertically grown mask
    // --------------------------------------------------------------------
    // Initialize grow volume
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(perimeter_data + i) == 1) {
            *(voronoi_data + i) = 1;
            *(flood_step_data + i) = 1.;
            *(flood_dist_data + i) = 1.;
        } else if (*(nii_midgm_data + i) == 1) {
            *(voronoi_data + i) = 2;
            *(flood_step_data + i) = 1.;
            *(flood_dist_data + i) = 1.;
        } else {
            *(voronoi_data + i) = 0;
            *(flood_step_data + i) = 0.;
            *(flood_dist_data + i) = 0.;
        }
    }

    grow_step = 1;
    voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t iii = 0; iii != nr_voi2; ++iii) {
            i = *(voi_id2 + iii);
            if (*(flood_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                bool jump_lock = false;
                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(voronoi_data + j) = *(voronoi_data + i);
                        }
                    } else if (*(nii_rim_data + j) != 0) {
                        jump_lock = true;
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (jump_lock == false) {

                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 3-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }

    // This is useful as a vertically grown rim mask
    for (uint32_t iii = 0; iii != nr_voi2; ++iii) {
        i = *(voi_id2 + iii);
        if (*(voronoi_data + i) == 2) {
            *(voronoi_data + i) = 0;
        }
    }
    save_output_nifti(fout, "perimeter_chunk", voronoi, false);

    // ========================================================================
    // Find quartets and compute distances on midgm domain
    // ========================================================================
    // Find first point on perimeter
    int idx_point1;
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);  // Map subset to full set
        if (*(perimeter_data + i) == 2) {
            idx_point1 = i;
        }
    }
    *(quartet_data + idx_point1) = 1;

    // Loop until desired number of points reached
    for (int32_t n = 2; n < 5; ++n) {
        int32_t grow_step = 1;
        uint32_t voxel_counter = nr_voxels;
        uint32_t ix, iy, iz, i, j;
        float d;

        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(quartet_data + i) != 0) {
                *(flood_step_data + i) = 1.;
                *(flood_dist_data + i) = 0.;
            } else {
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }

        // Reset some parameters
        grow_step = 1;
        voxel_counter = nr_voxels;
        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);  // Map subset to full set
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(perimeter_data + j) == 2) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                }
            }
            grow_step += 1;
        }

        // Find farthest point
        float max_distance = 0;
        int idx_new_point;
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);
            if (*(perimeter_data + i) == 2) {
                if (*(flood_dist_data + i) > max_distance) {
                    max_distance = *(flood_dist_data + i);
                    idx_new_point = i;
                }
            }
        }
        *(quartet_data + idx_new_point) = n;
    }

    if (mode_debug) {
        save_output_nifti(fout, "quartet", quartet, false);
    }

    // ========================================================================
    // Flooding distance from each quartet as output
    // ========================================================================
    for (int p = 1; p < 5; ++p) {
        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(quartet_data + i) == p) {
                *(flood_step_data + idx_point1) = 1.;
                *(flood_dist_data + idx_point1) = 1.;
            } else {
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }

        // Reset some parameters
        flood_dist_thr = std::numeric_limits<float>::infinity();
        grow_step = 1;
        voxel_counter = nr_voxels;

        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);  // Map subset to full set
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                }
            }
            grow_step += 1;
        }

        if (mode_debug) {
            save_output_nifti(fout, "point" + std::to_string(p) + "_midgm_dist", flood_dist, false);
        }

        // --------------------------------------------------------------------
        // Final Voronoi for propagating distances to all gray matter
        // --------------------------------------------------------------------
        cout << "\n  Start final Voronoi propagation..." << endl;
        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_midgm_data + i) != 0) {
                *(voronoi_data + i) = *(flood_dist_data + i);
                *(flood_step_data + i) = 1.;
                *(flood_dist_data + i) = 1.;
            } else {
                *(voronoi_data + i) = 0;
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }

        grow_step = 1;
        voxel_counter = nr_voxels;
        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t iii = 0; iii != nr_voi2; ++iii) {
                i = *(voi_id2 + iii);
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    bool jump_lock = false;
                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(nii_rim_data + j) == 3) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(voronoi_data + j) = *(voronoi_data + i);
                            }
                        } else if (*(nii_rim_data + j) != 0) {
                            jump_lock = true;
                        }
                    }

                    // ------------------------------------------------------------
                    // 2-jump neighbours
                    // ------------------------------------------------------------
                    if (jump_lock == false) {

                        if (ix > 0 && iy > 0) {
                            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xy;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iy < end_y) {
                            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xy;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy > 0) {
                            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xy;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy < end_y) {
                            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xy;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_yz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_yz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_yz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_yz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iz > 0) {
                            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iz > 0) {
                            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }

                        // ------------------------------------------------------------
                        // 3-jump neighbours
                        // ------------------------------------------------------------
                        if (ix > 0 && iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy > 0 && iz > 0) {
                            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix > 0 && iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy > 0 && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy < end_y && iz > 0) {
                            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                        if (ix < end_x && iy < end_y && iz < end_z) {
                            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                            if (*(nii_rim_data + j) == 3) {
                                d = *(flood_dist_data + i) + dia_xyz;
                                if (d < *(flood_dist_data + j)
                                    || *(flood_dist_data + j) == 0) {
                                    *(flood_dist_data + j) = d;
                                    *(flood_step_data + j) = grow_step + 1;
                                    *(voronoi_data + j) = *(voronoi_data + i);
                                }
                            }
                        }
                    }
                }
            }
            grow_step += 1;
        }
        // ----------------------------------------------------------------
        // Smooth distances
        // ----------------------------------------------------------------
        cout << "\n    Smoothing transitions..." << endl;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(smooth_data + i) = 0;
        }

        // Pre-compute weights
        float FWHM_val = 1;  // TODO(Faruk): Might tweak this one
        float w_0 = gaus(0, FWHM_val);
        float w_dX = gaus(dX, FWHM_val);
        float w_dY = gaus(dY, FWHM_val);
        float w_dZ = gaus(dZ, FWHM_val);

        for (uint16_t n = 0; n != 3; ++n) {
            for (uint32_t iii = 0; iii != nr_voi2; ++iii) {
                i = *(voi_id2 + iii);
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                float new_val = 0, total_weight = 0;

                // Start with the voxel itself
                new_val += *(voronoi_data + i) * w_0;
                total_weight += w_0;

                // --------------------------------------------------------
                // 1-jump neighbours
                // --------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dX;
                        total_weight += w_dX;
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dX;
                        total_weight += w_dX;
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dY;
                        total_weight += w_dY;
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dY;
                        total_weight += w_dY;
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dZ;
                        total_weight += w_dZ;
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        new_val += *(voronoi_data + j) * w_dZ;
                        total_weight += w_dZ;
                    }
                }
                *(smooth_data + i) = new_val / total_weight;
            }
            // Swap image data
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                *(voronoi_data + i) = *(smooth_data + i);
            }
        }

    save_output_nifti(fout, "point" + std::to_string(p) + "_dist", voronoi, true);
    }
    cout << "\n  Finished." << endl;
    return 0;
}
