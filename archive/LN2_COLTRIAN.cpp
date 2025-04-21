
#include "../dep/laynii_lib.h"
#include <fstream>

int shoq_help(void) {
    printf(
    "LN2_COLTRIAN: WORK IN PROGRESS... EXTREMELY EXPERIMENTAL!\n"
    "              Column triangulation.\n"
    "\n"
    "Usage:\n"
    "    LN2_COLTRIAN -input1 columns.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input1 : Mid-gray matter columns image from LN2_COLUMNS.\n"
    "    -input2 : Centroids file from LN2_COLUMNS.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    uint16_t ac;
    bool mode_debug = false;

    // Process user options
    if (argc < 2) return shoq_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return shoq_help();
        } else if (!strcmp(argv[ac], "-input1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input1\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-input2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input2\n");
                return 1;
            }
            fin2 = argv[ac];
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
        fprintf(stderr, "** missing option '-input1'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-input2'\n");
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

    log_welcome("LN2_COLTRIAN");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // ------------------------------------------------------------------------
    // Parse input file to determine output path
    string ext, basename;
    string path_out = fout;
    auto const pos = path_out.find_first_of('.');
    if (pos != string::npos) {
        basename = path_out.substr(0, pos);
    } else {  // Determine default extension when no extension given
        basename = path_out;
    }
    ext = "column_triangulated.obj";
    path_out = basename + "_" + ext;
    cout << "  Output:\n    " << path_out << endl;

    // Prepare to write wavefront .obj file
    ofstream myfile;
    myfile.open(path_out);

    // ------------------------------------------------------------------------
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

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_columns = copy_nifti_as_int32(nii1);
    int32_t* nii_columns_data = static_cast<int32_t*>(nii_columns->data);

    // ------------------------------------------------------------------------
    // Voxels of interest
    uint32_t nr_voi = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_columns_data + i) > 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_columns_data + i) > 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ------------------------------------------------------------------------
    // Prepare centroid coordinate array
    nifti_image* nii_centroids = copy_nifti_as_int32(nii2);
    int32_t* nii_centroids_data = static_cast<int32_t*>(nii_centroids->data);

    uint32_t nr_vertex = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_centroids_data + i) > 0){
            nr_vertex += 1;
        }
    }

    // Allocate 3D vertex coordinates
    int32_t* voi_coord;
    voi_coord = (int32_t*) malloc(nr_vertex*3*sizeof(int32_t));

    // Fill in coordinates
    // NOTE(Faruk): Strictly written for LN2_COLUMNS output columns. There
    // should be no integer skips for column IDs.
    uint32_t ix, iy, iz, j;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_centroids_data + i) > 0){
            int32_t column_id = *(nii_centroids_data + i) - 1;
            // Get array indices
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
            // Adjust to make voxel indices voxel centers
            float vx = static_cast<float>(ix) + 0.5;
            float vy = static_cast<float>(iy) + 0.5;
            float vz = static_cast<float>(iz) + 0.5;
            // Store float voxel center (or vertices)
            *(voi_coord + column_id * 3 + 0) = vx;
            *(voi_coord + column_id * 3 + 1) = vy;
            *(voi_coord + column_id * 3 + 2) = vz;
        }
    }

    // Write the vertices to file in order.
    for (uint32_t i = 0; i != nr_vertex; ++i) {
        myfile << "v " << *(voi_coord + i * 3 + 0) << " "
            << *(voi_coord + i * 3 + 1) << " "
            << *(voi_coord + i * 3 + 2) << "\n";
    }

    // ------------------------------------------------------------------------
    // Prepare triplets array to store potential triangular faces
    int32_t size_triplets = nr_voi * 3;
    int32_t* triplets;
    triplets = (int32_t*) malloc(size_triplets*sizeof(int32_t));
    for (int32_t i = 0; i != size_triplets; ++i) {
        *(triplets + i) = 0;
    }
    int count_triplets = 0;

    // ========================================================================
    // Triangulation on subset Voronoi cells
    // ========================================================================
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        // Keep neighbouring column ID's in an array
        int32_t* neighbours;
        neighbours = (int32_t*) malloc(8*sizeof(int32_t));
        for (int jj=0; jj<8; ++jj) {
            *(neighbours + jj) = 0;
        }
        *(neighbours + 0) = *(nii_columns_data + i);

        // 1-jump neighbours
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 1) = *(nii_columns_data + j);
            }
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 2) = *(nii_columns_data + j);
            }
        }
        if (iz > 0) {
            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 3) = *(nii_columns_data + j);
            }
        }

        // 2-jump neighbours
        if (ix > 0 && iy > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 4) = *(nii_columns_data + j);
            }
        }
        if (iy > 0 && iz > 0) {
            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 5) = *(nii_columns_data + j);
            }
        }
        if (ix > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 6) = *(nii_columns_data + j);
            }
        }

        // 3-jump neighbour
        if (ix > 0 && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            if (*(nii_columns_data + j) > 0) {
                *(neighbours + 7) = *(nii_columns_data + j);
            }
        }

        // --------------------------------------------------------------------
        // Bubble sort
        int m, n, flag = 1;  // set flag to 1 to start first pass
        int temp;  // holding variable
        for(m = 1; (m <= 8) && flag; m++) {
            flag = 0;
            for (n=0; n < (8 - 1); n++) {
                if (*(neighbours + n + 1) < *(neighbours + n)) {
                    temp = *(neighbours + n);  // swap elements
                    *(neighbours + n) = *(neighbours + n + 1);
                    *(neighbours + n + 1) = temp;
                    flag = 1; // indicates that a swap occurred.
                }
            }
        }

        // --------------------------------------------------------------------
        // Find non-zero unique elements in array
        int count_neigh = 0;

        // Have an array to note down columns id's that will form triangles.
        // NOTE(Faruk): I allocate 8 elements here because if the volume is
        // over-parcellated (i.e. too many columns), hit count can go above 3
        int32_t* trio;
        trio = (int32_t*) malloc(8*sizeof(int32_t));

        for (int m=0; m<8; m++) {
            int n;
            for (n=0; n<m; n++) {
                if (*(neighbours + m) == *(neighbours + n)) {
                    break;
                }
            }
            if (m == n) {
                if (*(neighbours + m) != 0) {
                    *(trio + count_neigh) = *(neighbours + m);
                    count_neigh += 1;
                }
            }
        }
        if (count_neigh == 3) {  // Hard-lock to triplets only
            for (int l=0; l<count_neigh; ++l) {
                *(triplets + count_triplets * 3 + l) = *(trio + l);
            }
            count_triplets += 1;
        }
    }

    // Write triples (faces) into wavefront .obj file
    for (int i=0; i<count_triplets; ++i) {
        myfile << "f ";
        for (int j=0; j<3; ++j){
            myfile << *(triplets + i * 3 + j) << " ";
        }
        myfile << "\n";
    }

    myfile.close();

    // save_output_nifti(fout, "gradients", nii_grad, true);

    cout << "\n  Finished." << endl;
    return 0;
}
