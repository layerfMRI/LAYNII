#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN2_COPY_GEOMETRY: Copies the spatial geometry (voxel size, qform, sform\n"
    "                   and spatial units) from the header of a reference\n"
    "                   nifti file to another nifti file. Only the local\n"
    "                   geometry is copied. All other header information of the\n"
    "                   input file (data type, number of time points, TR,\n"
    "                   scaling, description, etc.) and its data are kept.\n"
    "\n"
    "Usage:\n"
    "    LN2_COPY_GEOMETRY -input destination.nii -geometry origin.nii\n"
    "    ../LN2_COPY_GEOMETRY -input lo_BOLD_intemp.nii -geometry lo_Nulled_intemp.nii\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -input    : Nifti file whose header geometry will be replaced.\n"
    "                Its data and non-geometry header fields are kept.\n"
    "    -geometry : Reference nifti file from which the geometry is taken.\n"
    "                Must have the same number of voxels along x, y, and z\n"
    "                as the input. It can have a different data type and a\n"
    "                different number of time points.\n"
    "    -output   : (Optional) Output filename, including .nii or\n"
    "                .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n"
    "Notes:\n"
    "    The following fields are copied from the reference:\n"
    "        pixdim[1-3] (dx, dy, dz), qform_code, quatern_b/c/d,\n"
    "        qoffset_x/y/z, qfac (pixdim[0]), sform_code, srow_x/y/z,\n"
    "        and the spatial part of xyzt_units.\n"
    "    The following fields are NOT copied (taken from the input):\n"
    "        dim (incl. number of time points), datatype, bitpix,\n"
    "        pixdim[4-7] (incl. TR), time units, toffset, scl_slope,\n"
    "        scl_inter, cal_min, cal_max, dim_info, intent, descrip.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL, *fgeom = NULL, *fout = NULL;
    bool use_outpath = false;
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
            fin = argv[ac];
            if (!use_outpath) fout = argv[ac];
        } else if (!strcmp(argv[ac], "-geometry")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -geometry\n");
                return 1;
            }
            fgeom = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
            use_outpath = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fgeom) {
        fprintf(stderr, "** missing option '-geometry'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }
    // Read reference dataset, header only (data is not needed)
    nifti_image* nii_geom = nifti_image_read(fgeom, 0);
    if (!nii_geom) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fgeom);
        return 2;
    }

    log_welcome("LN2_COPY_GEOMETRY");
    log_nifti_descriptives(nii_input);
    log_nifti_descriptives(nii_geom);

    // ========================================================================
    // Sanity check: the spatial matrix size must match
    // ========================================================================
    if (nii_input->nx != nii_geom->nx || nii_input->ny != nii_geom->ny
        || nii_input->nz != nii_geom->nz) {
        fprintf(stderr, "** spatial dimensions do not match:\n"
                "   input    : %lld x %lld x %lld\n"
                "   geometry : %lld x %lld x %lld\n",
                (long long)nii_input->nx, (long long)nii_input->ny,
                (long long)nii_input->nz, (long long)nii_geom->nx,
                (long long)nii_geom->ny, (long long)nii_geom->nz);
        return 3;
    }

    // ========================================================================
    // Copy only the local (spatial) geometry
    // ========================================================================
    cout << "  Copying geometry..." << endl;

    // Voxel sizes
    nii_input->dx = nii_geom->dx;
    nii_input->dy = nii_geom->dy;
    nii_input->dz = nii_geom->dz;
    nii_input->pixdim[1] = nii_geom->pixdim[1];
    nii_input->pixdim[2] = nii_geom->pixdim[2];
    nii_input->pixdim[3] = nii_geom->pixdim[3];

    // Spatial units only, time units of the input are kept
    nii_input->xyz_units = nii_geom->xyz_units;

    // Qform (quaternion representation and derived matrices)
    nii_input->qform_code = nii_geom->qform_code;
    nii_input->quatern_b = nii_geom->quatern_b;
    nii_input->quatern_c = nii_geom->quatern_c;
    nii_input->quatern_d = nii_geom->quatern_d;
    nii_input->qoffset_x = nii_geom->qoffset_x;
    nii_input->qoffset_y = nii_geom->qoffset_y;
    nii_input->qoffset_z = nii_geom->qoffset_z;
    nii_input->qfac = nii_geom->qfac;
    nii_input->pixdim[0] = nii_geom->pixdim[0];
    nii_input->qto_xyz = nii_geom->qto_xyz;
    nii_input->qto_ijk = nii_geom->qto_ijk;

    // Sform (affine matrix and its inverse)
    nii_input->sform_code = nii_geom->sform_code;
    nii_input->sto_xyz = nii_geom->sto_xyz;
    nii_input->sto_ijk = nii_geom->sto_ijk;

    // ========================================================================
    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "geomcopied", nii_input, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
