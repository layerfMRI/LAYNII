
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "../dep/nifti2_io.h"
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_CONLAY: This program consenses the layers on a lower resolution grid \n"
    "           This can be usefull, when the layerification has been  .\n"
    "           originally done on an upscaled anatomical image and you \n"
    "           want to use it on lower resolution functional data \n"
    "           This program has been originally written for Federico \n"
    "\n"
    "           The two used spatial grids need to be scaled by intager values \n"
    "           this program does not have a interpolation adapter \n"
    "\n"
    "Usage: \n"
    "    LN_CONLAY -layers highres_layers.nii -ref funct.nii -output lowres_layers.nii\n"
    "    ../LN_CONLAY -layers lo_sc_layers.nii -ref lo_T1EPI.nii -subsample -output lo_layers_out.nii.gz \n"
    "\n"
    "Options:\n"
    "\n"
    "    -help       : Show this help.\n"
    "    -layers     : Nifti (.nii) file that should be zoomed (e.g. with \n"
    "                     multiple time points).\n"
    "    -ref        : Nifti (.nii) file that determines the region of interest\n"
    "                     (e.g. the layer mask with one time point).\n"
    "    -output     : (Optional) Output filename, including .nii or\n"
    "                  .nii.gz, and path if needed. Overwrites existing files.\n"    
    "    -subsample  : This option is regridds the layer values based on the voxels centroid\n"
    "                  If this option is not used (default, the program \n"
    "                  takes the local average instead \n"
    "                  I would only recommend this for odd scale factors \n"
    "                  otherwiese your results will be dependent on the \n"
    "                  specific convention of upscaling tools e.g. AFNI!=BV \n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    bool use_outpath = false ;
    bool subsample = false ;
    char *fout = NULL ;
    char *fin_1 = NULL, *fin_2 = NULL;
    int ac;
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-ref")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -ref\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-subsample")) {
            subsample = true;
            cout << "I am subsampling the voxel centroid"  << endl;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-ref'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_CONLAY");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);  // ref

    // Get dimensions of input
    const int size_z = nii1->nz;
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_time = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // get high res layers
    nifti_image* nim_file_1 = copy_nifti_as_int16(nii1);
    short* nii1_data = static_cast<short*>(nim_file_1->data);

    // ========================================================================
    // Get grid size for low res layers

    //nifti_image* nii_outlay = copy_nifti_as_float32(nii2);
    //float* nii_outlay_data = static_cast<float*>(nii_outlay->data);

    nifti_image* nii_outlay = nifti_copy_nim_info(nii2);
    nii_outlay->datatype = NIFTI_TYPE_FLOAT32;
    nii_outlay->nbyper = sizeof(float);
    nii_outlay->nvox =  nii_outlay->nvox/nii_outlay->nt ;
    nii_outlay->nt = 1;  // layers only have one time point
    nii_outlay->data = calloc(nii_outlay->nvox, nii_outlay->nbyper);
    float* nii_outlay_data = static_cast<float*>(nii_outlay->data);
    int nrout_voxels = nii_outlay->nvox;

    int size_zout = nii_outlay->nz;
    int size_xout = nii_outlay->nx;
    int size_yout = nii_outlay->ny;
    int size_timeout = nii_outlay->nt;
    const int nxo = nii_outlay->nx;
    const int nxyo = nii_outlay->nx * nii_outlay->ny;
    const int nxyzo = nii_outlay->nx * nii_outlay->ny * nii_outlay->nz;
    nii_outlay->scl_slope = nii1->scl_slope ; // To make sure that the layers are in the same units


    // ========================================================================
    // Get ratios of resolutions and matrix sizes
    double sratio_x = (double)size_x/(double)size_xout ;
    double sratio_y = (double)size_y/(double)size_yout ;
    double sratio_z = (double)size_z/(double)size_zout ;
    cout << "    Matrix size reduction  =  |" << sratio_x << "X|  x  |" << sratio_y   << "Y| x |" << sratio_z << "Z| " <<  endl;

    double rratio_x = (double)(nii_outlay->pixdim[1])/(double)(nii1->pixdim[1]) ;
    double rratio_y = (float)(nii_outlay->pixdim[2])/(double)(nii1->pixdim[2]) ;
    double rratio_z = (float)(nii_outlay->pixdim[3])/(double)(nii1->pixdim[3]) ;

    cout << "    Voxel size reduction   =  |" << rratio_x << "X|  x  |" << rratio_y   << "Y| x |" << rratio_z << "Z| " << endl;

    // ========================================================================
    // complaining if something is weird.
    if ((fmod(sratio_x,1)>4.76838e-06 && fmod(sratio_x,1)<1.-4.76838e-06)|| (fmod(sratio_y,1)>4.76838e-06 && fmod(sratio_y,1)<1.-4.76838e-06 ) || (fmod(sratio_z,1)>4.76838e-06 && fmod(sratio_z,1)<1.-4.76838e-06)) { //4.76838e-07is precision of float
        cout << " ******************************************* " << endl;
        cout << " *** there is a non intager ratio of the *** " << endl;
        cout << " *** matrix size. This is not recommended ** " << endl;
        cout << " *** I hope you know what you are doing. *** " << endl;
        cout << " ******************************************* " << endl;
    }
    if ((fmod(rratio_x,1)>4.76838e-06 && fmod(rratio_x,1)<1.-4.76838e-06)|| (fmod(rratio_y,1)>4.76838e-06 && fmod(rratio_y,1)<1.-4.76838e-06 ) || (fmod(rratio_z,1)>4.76838e-06 && fmod(rratio_z,1)<1.-4.76838e-06)) {
        cout << " ******************************************* " << endl;
        cout << " *** there is a non intager ratio of the *** " << endl;
        cout << " *** voxel size. This is not recommended *** " << endl;
        cout << " *** I hope you know what you are doing. *** " << endl;
        cout << " ******************************************* " << endl;
    }

    if ( (rratio_x - sratio_x)>4.76838e-06 || (rratio_y - sratio_y)>4.76838e-06 || (rratio_z - sratio_z)>4.76838e-06 || (rratio_x - sratio_x)<-4.76838e-06 || (rratio_y - sratio_y)<-4.76838e-06 || (rratio_z - sratio_z)<-4.76838e-06 ) {
        cout << " ******************************************* " << endl;
        cout << " *** The matrix size and the voxel side  *** " << endl;
        cout << " *** do not match across datasets.       *** " << endl;
        cout << " *** If your headers are right, these    *** " << endl;
        cout << " *** cannot be used here.                *** " << endl;
        cout << " *** If your headers are not correct     *** " << endl;
        cout << " *** Please fix them first.              *** " << endl;
        cout << " ******************************************* " << endl;
    }



    // ========================================================================
    // Copy voxel from the bigger input nifti image by averaging
    if (!subsample) {
        cout << "    I am averaging not subsampling" << endl;
        int lo_z, lo_y, lo_x ;
        float mean_val = 0.;
        float nr_vox = 0;


        for (int it = 0; it < size_timeout; ++it) {
            for (int iz = 0; iz <size_z ; iz=iz+sratio_z) {
                for (int iy = 0 ; iy < size_y; iy=iy+sratio_y) {
                    for (int ix = 0; ix < size_x; ix=ix+sratio_x) {
                        mean_val = 0.;
                        nr_vox = 0;
                        for (int jz = iz; jz < iz+sratio_z ; ++jz) {
                            for (int jy = iy ; jy < iy+sratio_y; ++jy) {
                                for (int jx = ix; jx < ix+sratio_x; ++jx) {
                                    if ( *(nii1_data + nxyz * it + nxy * jz + nx * jy + jx) > 0 ){
                                        mean_val = mean_val + (float)(*(nii1_data + nxyz * it + nxy * jz + nx * jy + jx)) ;
                                        nr_vox++;
                                    }
                                }
                            }
                        }

                        lo_x = ix/sratio_x ;
                        lo_y = iy/sratio_y ;
                        lo_z = iz/sratio_z ;
                        if (nr_vox != 0 ) mean_val = mean_val/nr_vox ;
                        *(nii_outlay_data + nxyzo * it + nxyo * lo_z + nxo * lo_y + lo_x) =  mean_val;
                    }
                }
            }
        }
    }


    // ========================================================================
    // Copy voxel from the bigger input nifti image by subsampling

    if (subsample) {
        cout << "I am subsampling not averaging " << endl;
        int up_z, up_y, up_x ;

        for (int it = 0; it < size_timeout; ++it) {
            for (int iz = 0; iz <size_zout ; ++iz) {
                for (int iy = 0 ; iy < size_yout; ++iy) {
                    for (int ix = 0; ix < size_xout; ++ix) {
                        up_x = ix*sratio_x + (int)(sratio_x)/2 ; // The division with two is to calculate based on voxel centroits rather than voxel edges
                        up_y = iy*sratio_y + (int)(sratio_y)/2 ;
                        up_z = iz*sratio_z + (int)(sratio_z)/2 ;
                        *(nii_outlay_data + nxyzo * it + nxyo * iz + nxo * iy + ix) =  (float)(*(nii1_data + nxyz * it + nxy * up_z + nx * up_y + up_x)) ;
                    }
                }
            }
        }
    }

    if (!use_outpath) fout = fin_2;
    save_output_nifti(fout, "sub_layers", nii_outlay, true, use_outpath);

    cout << "Finished!" << endl;
    return 0;
}
