#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_HEXBIN: Generate hexagonal bins based on UV coordinates generated \n"
    "            from LN2_MULTILATERATE.\n"
    "\n"
    "Usage:\n"
    "    LN2_HEXBIN -coord_uv coord_uv.nii -radius 10\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -coord_uv : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -radius   : Radius of the circle inscribed within hexagons.\n"
    "                In UV coordinate metric units (e.g. mm)."
    "    -output   : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    float radius = 10;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-coord_uv")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_uv\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-radius")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radius\n");
                return 1;
            }
            radius = atof(argv[ac]);
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
        fprintf(stderr, "** missing option '-coords_uv'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_HEXBIN");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const int nr_voxels = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // nifti_image* nii_bins = copy_nifti_as_int32(nii1);
    // float* nii_bins_data = static_cast<float*>(nii_bins->data);

    // ========================================================================
    // Allocating new nifti for output 3D nifti
    nifti_image* nii_bins = nifti_copy_nim_info(nii_input);
    nii_bins->datatype = NIFTI_TYPE_INT32;
    nii_bins->dim[0] = 4;  // For proper 4D nifti
    // nii_bins->dim[1] = bins_u;
    // nii_bins->dim[2] = bins_v;
    // nii_bins->dim[3] = bins_d;
    nii_bins->dim[4] = 1;
    nifti_update_dims_from_array(nii_bins);
    nii_bins->nvox = nr_voxels;
    nii_bins->nbyper = sizeof(int32_t);
    nii_bins->data = calloc(nii_bins->nvox, nii_bins->nbyper);
    nii_bins->scl_slope = 1;
    int32_t* nii_bins_data = static_cast<int32_t*>(nii_bins->data);

    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_bins_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    int nr_voi = 0;  // Voxels of interest
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int* voi_id;
    voi_id = (int*) malloc(nr_voi*sizeof(int));

    // Fill in indices to be able to remap from subset to full set of voxels
    int ii = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find coordinate ranges
    // ========================================================================
    float min_u = std::numeric_limits<float>::max();
    float max_u = std::numeric_limits<float>::min();
    float min_v = std::numeric_limits<float>::max();
    float max_v = std::numeric_limits<float>::min();

    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        // Check U coordinate min & max
        if (*(nii_input_data + nr_voxels*0 + i) < min_u) {
            min_u = *(nii_input_data + nr_voxels*0 + i);
        }
        if (*(nii_input_data + nr_voxels*0 + i) > max_u) {
            max_u = *(nii_input_data + nr_voxels*0 + i);
        }
        // Check V coordinate min & max
        if (*(nii_input_data + nr_voxels*1 + i) < min_v) {
            min_v = *(nii_input_data + nr_voxels*1 + i);
        }
        if (*(nii_input_data + nr_voxels*1 + i) > max_v) {
            max_v = *(nii_input_data + nr_voxels*1 + i);
        }
    }
    cout << "  U coordinate min & max: " << min_u << " | " << max_u << endl;
    cout << "  V coordinate min & max: " << min_v << " | " << max_v << endl;

    // ========================================================================
    // Place hexbin centers within UV range
    // ========================================================================
    // Figure out orthogonal step sizes
    float diameter = radius * 2;
    float step_u = diameter;
    float step_v = diameter / sqrt(2);

    // Guesstimate number of bins needed
    int nr_bins_u = (abs(min_u) + abs(max_u)) / step_u;
    int nr_bins_v = (abs(min_v) + abs(max_v)) / step_v;
    int nr_bins = nr_bins_u * nr_bins_v;

    // Initialize all elements at 0
    float* arr_centers_u = new float[nr_bins]();
    float* arr_centers_v = new float[nr_bins]();

    // Set first center as minima
    arr_centers_u[0] = min_u;
    arr_centers_v[0] = min_v;

    // Compute hexbins centers
    int k;
    for (int j = 0; j != nr_bins_v; ++j) {
        for (int i = 0; i != nr_bins_u; ++i) {
            k = j * nr_bins_u + i;
            if (j % 2 == 0) {  // even rows
                arr_centers_u[k] = min_u + step_u * i;
                arr_centers_v[k] = min_v + step_v * j;
            } else {
                arr_centers_u[k] = min_u + (step_u * i) + (step_u / 2);
                arr_centers_v[k] = min_v + step_v * j;
            }
        }
    }

    // ========================================================================
    // Evaluate each voxel agains hexbin centers to find closest center
    // ========================================================================
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);

        float coord_u = *(nii_input_data + nr_voxels*0 + i);
        float coord_v = *(nii_input_data + nr_voxels*1 + i);

        // Compute distance to each bin center
        float min_dist = std::numeric_limits<float>::max();
        for (int32_t j = 0; j != nr_bins; j++) {
            float bin_u = arr_centers_u[j];
            float bin_v = arr_centers_v[j];

            float dist = sqrt(pow(coord_u - bin_u, 2) + pow(coord_v - bin_v, 2));
            if (dist < min_dist) {
                min_dist = dist;
                *(nii_bins_data + i) = j;
            }
        }
    }

    std::ostringstream tag;
    tag << radius;
    save_output_nifti(fout, "hexbins"+tag.str(), nii_bins, true);

    cout << "\n  Finished." << endl;
    return 0;
}
