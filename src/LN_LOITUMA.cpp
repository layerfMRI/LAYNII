#include "../dep/laynii_lib.h"

int show_help( void )
{
   printf(
      "LN_LOITUMA : This program generates equi-volume layers based on Leaky layers and equi-dist layers\n"
      "\n"
      "    This program does not necessarily assume any curcatures smoothness ,\n"
      "    ,\n"
      "\n"
      "    basic usage: LN_LOITUMA -equidist layers_equi_dist.nii -leaky layers_leaky.nii \n"
      "\n"
      "\n"
      "   This program now supports INT16, INT32 and FLOAT23 \n"

      "\n"
      "       -help                 : show this help\n"
      "       -equidist             : nii file that contains many many layers \n"
      "                               Ideally, the cortex is devided into 1000 layers.  \n"
      "                               The layers are estimated based on the equi-distance principle.  \n"
      "       -leaky                : nii file that contains many many layers \n"
      "                               Ideally, the cortex is devided into 1000 layers.  \n"
      "                               The layers are estimated based on the leaky-layer principle.  \n"
      "       -FWHM                 : Optional parameter to enforce a smooth curvature \n"
      "       -nr_layers            : Optional parameter to enforce a smooth curvature \n"
      "       -output               : Optional parameter for output fine name (path) \n"
      "                               default is equivol_layers_.nii in current folder \n"
      "                                \n"
      "                                \n"
      "                                \n"
      "                               If you run this on EPI-T1 data consider preparing them as follwos, E.g:  \n"
      "                               LN_GROWLAYERS -rim rim.nii -N 1000 \n"
      "                               LN_LEAKY_LAYERS -rim rim.nii -nr_layers 100 \n"     
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{
   bool use_outpath = false;
   char       * fleakyi=NULL, * fdisti=NULL,  * froii=NULL ;
   const char * fout="equivol_layers.nii";
   int          ac ; 
   float        FWHM_val=1  ;
   if( argc < 3 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-equidist") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -equidist\n");
            return 1;
         }
         fdisti = argv[ac];  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-FWHM") ) {
        if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -FWHM\n");
            return 1;
         }
         FWHM_val = atof(argv[ac]);  // no string copy, just pointer assignment 
         cout << "I will assume a smooth folding pattern"  << endl; 
         //FWHM_val = 1 ; // do smoothing anyway 
      }
      else if( ! strcmp(argv[ac], "-leaky") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -leaky\n");
            return 1;
         }
         fleakyi = argv[ac];  // no string copy, just pointer assignment 
      }
      else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
      } 
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !fleakyi  ) { fprintf(stderr, "** missing option '-leaky'\n");  return 1; }
//   // read input dataset, including data 
//   nifti_image * nim_leakyi = nifti_image_read(fleakyi, 1);
//   if( !nim_leakyi ) {
//      fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", fleakyi);
//      return 2;
//   }
   
   if( !fdisti  ) { fprintf(stderr, "** missing option '-equidist'\n");  return 1; }
   // read input dataset, including data 
//   nifti_image * nim_disti = nifti_image_read(fdisti, 1);
//   if( !nim_disti ) {
//      fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", fdisti);
//      return 2;
//   }
   
    nifti_image* fdist = nifti_image_read(fdisti, 1);
    nifti_image* nii_dist = copy_nifti_as_int16(fdist);
    int16_t* nii_dist_data = static_cast<int16_t*>(nii_dist->data);
   
    nifti_image* fleaky = nifti_image_read(fleakyi, 1);
    nifti_image* nii_leak = copy_nifti_as_int16(fleaky);
    int16_t* nii_leak_data = static_cast<int16_t*>(nii_leak->data);
    
    
     // Get dimensions of input
    const uint32_t size_z = nii_leak->nz;
    const uint32_t size_x = nii_leak->nx;
    const uint32_t size_y = nii_leak->ny;
    const uint32_t nrep =  nii_leak->nt; 
    const uint32_t nx = nii_leak->nx;
    const uint32_t nxy = nii_leak->nx * nii_leak->ny;
    const uint32_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii_leak->pixdim[1];
    const float dY = nii_leak->pixdim[2];
    const float dZ = nii_leak->pixdim[3];
    

    log_welcome("LN_LOITUMA");
    log_nifti_descriptives(nii_leak);
    log_nifti_descriptives(nii_dist);



   if (FWHM_val > 0 ){
       uint32_t iter_smooth = 1; 
       uint32_t ix, iy, iz, j, k;
       float w = 0;
       

        // ------------------------------------------------------------------------
        // Smooth equi-dist data within GM
        // ------------------------------------------------------------------------
        nifti_image* smooth = copy_nifti_as_float32(nii_leak);
        float* smooth_data = static_cast<float*>(smooth->data);
        nifti_image* smooth_eqdist = copy_nifti_as_float32(nii_leak);
        float* smooth_eqdist_data = static_cast<float*>(smooth_eqdist->data);
        
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(smooth_data + i) = 0;
            *(smooth_eqdist_data + i ) = *(nii_dist_data + i ) ; 
        }

        for (uint16_t n = 0; n != iter_smooth; ++n) {
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(nii_leak_data + i) > 0) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    float new_val = 0, total_weight = 0;

                    // Start with the voxel itself
                    w = gaus(0, FWHM_val);
                    new_val += *(smooth_eqdist_data + i) * w;
                    total_weight += w;

                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix != 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dX, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (ix != size_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dX, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iy != 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dY, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iy != size_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dY, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iz != 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dZ, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iz != size_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dZ, FWHM_val);
                            new_val += *(smooth_eqdist_data + j) * w;
                            total_weight += w;
                        }
                    }
                    *(smooth_data + i) = new_val / total_weight;
                }
            }
            // Swap image data
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                *(smooth_eqdist_data + i) = *(smooth_data + i);
            }
        }

            save_output_nifti(fout, "smooth_equi_dist.nii", smooth_eqdist, false);
        for (uint32_t i = 0; i != nr_voxels; ++i) *(nii_dist_data + i) = *(smooth_eqdist_data + i) ; 
   } // smoothing loop closed
   

        // ------------------------------------------------------------------------
        // estimating equ vol
        // ------------------------------------------------------------------------

        nifti_image* equi_vol = copy_nifti_as_int16(nii_leak);
        short* equi_vol_data = static_cast<short*>(equi_vol->data);
        
int max_layer = 0; 
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if ( *(nii_dist_data + i ) > max_layer) max_layer = *(nii_dist_data + i )  ; 
    }
int min_layer = 100000; 
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if ( *(nii_dist_data + i) < min_layer && *(nii_dist_data + i) > 0 ) min_layer = *(nii_dist_data + i )  ; 
    }
    
cout << " min layer is " <<  min_layer << "   max layers is " << max_layer << endl; 

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_dist_data + i) > 0) {
            *(equi_vol_data + i) = *(nii_dist_data + i)* 2. - *(nii_leak_data+ i); 
            if (*(equi_vol_data + i)  > max_layer) *(equi_vol_data + i) = max_layer ; 
            if (*(equi_vol_data + i)  < min_layer) *(equi_vol_data + i) = min_layer ; 
            if (*(nii_dist_data + i) == max_layer) *(equi_vol_data + i) = max_layer ; 
            if (*(nii_dist_data + i) == min_layer) *(equi_vol_data + i) = min_layer ; 

           // if ( *(nii_leak_data + i) <= 0  ) *(equi_vol_data + i) = 0 ; 
        }
        else {
            *(equi_vol_data + i) = 0 ; 
        }
    }


        // ------------------------------------------------------------------------
        // cleaning up equivol
        // ------------------------------------------------------------------------



   if (FWHM_val > 0 ){
       uint32_t iter_smooth = 1; 
       uint32_t ix, iy, iz, j, k;
       float w = 0;
       
        nifti_image* smooth = copy_nifti_as_float32(nii_leak);
        float* smooth_data = static_cast<float*>(smooth->data);
        nifti_image* smooth_eqvol = copy_nifti_as_float32(nii_leak);
        float* smooth_eqvol_data = static_cast<float*>(smooth_eqvol->data);

        // ------------------------------------------------------------------------
        // Smooth equi-dist data within GM
        // ------------------------------------------------------------------------

        
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(smooth_data + i) = 0;
            *(smooth_eqvol_data + i ) = *(equi_vol_data + i ) ; 
        }

        for (uint16_t n = 0; n != iter_smooth; ++n   ) {
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(nii_leak_data + i) > 0) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    float new_val = 0, total_weight = 0;

                    // Start with the voxel itself
                    w = gaus(0, FWHM_val);
                    new_val += *(smooth_eqvol_data + i) * w;
                    total_weight += w;

                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix != 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dX, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (ix != size_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dX, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iy != 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dY, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iy != size_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dY, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iz != 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dZ, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    if (iz != size_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_leak_data + j) > 0) {
                            w = gaus(dZ, FWHM_val);
                            new_val += *(smooth_eqvol_data + j) * w;
                            total_weight += w;
                        }
                    }
                    *(smooth_data + i) = new_val / total_weight;
                }
            }
            // Swap image data
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                *(smooth_eqvol_data + i) = *(smooth_data + i);
            }
        }

            save_output_nifti(fout, "", smooth_eqvol,true, use_outpath);


   } // smoothing loop closed





    save_output_nifti(fout, "unsmoothed_equi_vol.nii", equi_vol,false); 

  return 0;
}


