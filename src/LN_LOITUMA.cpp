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
   float        FWHM_val=0  ;
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
/*
   
// if you are running the smoothing in 2D, it will still go thought he entire pipeline. 
// the only difference is that the weights in a certain direction are suppressed
// doing it in 2Dim, will not speed up the program
   if  (twodim == 1) dZ = 1000 * dZ ; 
   
   nifti_image * nim_inputf  	= nifti_copy_nim_info(nim_inputfi);
   nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
   nim_inputf->nbyper = sizeof(float);
   nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
   float  *nim_inputf_data = (float *) nim_inputf->data;
   
   nifti_image * nim_mask  	= nifti_copy_nim_info(nim_maski);
   nim_mask->datatype = NIFTI_TYPE_FLOAT32;
   nim_mask->nbyper = sizeof(float);
   nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
   float  *nim_mask_data = (float *) nim_mask->data;
   
   nifti_image * nim_roi  	= nifti_copy_nim_info(nim_maski);
   nim_roi->datatype = NIFTI_TYPE_FLOAT32;
   nim_roi->nbyper = sizeof(float);
   nim_roi->data = calloc(nim_roi->nvox, nim_roi->nbyper);
   float  *nim_roi_data = (float *) nim_roi->data; 
 
   /////////////////////////////////////////////////////////////////////////
   /////////  fixing potential problems with different input datatypes /////
   /////////  here, I am loading them in their native datatype /////////////
   /////////  and translate them to the datatime I like best  //////////////
   /////////////////////////////////////////////////////////////////////////

if ( nim_inputfi->datatype == NIFTI_TYPE_FLOAT32 ||  nim_inputfi->datatype == NIFTI_TYPE_INT32 ) {
  float  *nim_inputfi_data = (float *) nim_inputfi->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_inputf_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_inputfi_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}  
  

if ( nim_inputfi->datatype == NIFTI_TYPE_INT16 || nim_inputfi->datatype == DT_UINT16 ) {
  short  *nim_inputfi_data = (short *) nim_inputfi->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_inputf_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_inputfi_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    

if ( nim_inputfi->datatype == DT_FLOAT64 || nim_inputfi->datatype == NIFTI_TYPE_FLOAT64 ) {
  double  *nim_inputfi_data = (double *) nim_inputfi->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_inputf_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_inputfi_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    


if ( nim_maski->datatype == NIFTI_TYPE_FLOAT32 || nim_maski->datatype == NIFTI_TYPE_INT32  ) {
  float  *nim_maski_data = (float *) nim_maski->data;
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data   +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_maski_data  +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
}    
  
if ( nim_maski->datatype == NIFTI_TYPE_INT16 || nim_maski->datatype ==  DT_UINT16) {
  short  *nim_maski_data = (short *) nim_maski->data;
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data   +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_maski_data   +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
}    

if ( nim_maski->datatype == DT_FLOAT64 || nim_maski->datatype ==  NIFTI_TYPE_FLOAT64) {
  double  *nim_maski_data = (double *) nim_maski->data;
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data   +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_maski_data   +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
}    


if ( do_masking == 1 ) {
	
	
            // read input dataset, including data 
         nifti_image * nim_roii = nifti_image_read(froii, 1);
         if( !nim_roii ) {
           fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", froii);
           return 2;
         }	

 if ( nim_roii->datatype == NIFTI_TYPE_FLOAT32 ||  nim_roii->datatype ==  NIFTI_TYPE_INT32 ) {
  float  *nim_roii_data = (float *) nim_roii->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_roi_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_roii_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
 }  

 if ( nim_roii->datatype == NIFTI_TYPE_INT16 || nim_roii->datatype == DT_UINT16 ) {
  short  *nim_roii_data = (short *) nim_roii->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_roi_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_roii_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
 }    
}
	// write out some stuff that might be good to know, if you want to debug
   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 
   cout << " Voxel size    " <<  dX << " x " <<  dY << " x "  <<  dZ  << endl; 

	cout << " datatye 1 = " << nim_inputf->datatype << endl;
    cout << " datatye 2 = " << nim_mask ->datatype << endl;

/////////////////////////////////////////////
////   MAKE allocating necessary files  /////
/////////////////////////////////////////////

    
    nifti_image * smoothed  	= nifti_copy_nim_info(nim_inputf);
    nifti_image * gausweight  	= nifti_copy_nim_info(nim_inputf);
    smoothed->datatype 		= NIFTI_TYPE_FLOAT32; 
	gausweight->datatype 	= NIFTI_TYPE_FLOAT32;
    smoothed->nbyper 		= sizeof(float);
	gausweight->nbyper 		= sizeof(float);
    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
    
    // the gaus wieght is a gemoetry factor and only needs to be estimated once (not for every time step)
    // so with the next lines I am saving RAM
    gausweight->nt 				= 1	; 
    gausweight->nvox 			= gausweight->nvox / nrep ; 
    float  *smoothed_data = (float *) smoothed->data;
    float  *gausweight_data = (float *) gausweight->data;    

//    nifti_image * debug  	= nifti_copy_nim_info(nim_inputf);
//    debug->datatype 		= NIFTI_TYPE_FLOAT32; 
//    debug->nbyper 			= sizeof(float);
//    debug->data 			= calloc(debug->nvox, debug->nbyper);
//    float  *debug_data 		= (float *) debug->data;
    
 //   if ( do_masking == 1 ) {
//	  for(int islice=0; islice<sizeSlice; ++islice){  
//	      for(int iy=0; iy<sizePhase; ++iy){
//	        for(int ix=0; ix<sizeRead; ++ix){
  //      		 	*(debug_data +  nxy*islice + nx*ix  + iy  ) = 0; 
//           } 
//	    }
//	  }
//    }   


cout << " time dimension of soothed. output file:  " << smoothed->nt <<  endl; 



//float kernal_size = 10; // corresponds to one voxel sice. 
int vinc = max(1.,2. * FWHM_val/dX ); // if voxel is too far away, I ignore it. 
float dist_i = 0.;
cout << " vinc " <<  vinc<<  endl; 
cout << " FWHM_val " <<  FWHM_val<<  endl; 

float temp_wight_factor = 0 ; // to store temp values, so I don't make the same computations over and over again. 

///////////////////////////////////////////////////
////finding the range of gradient values   ////////
///////////////////////////////////////////////////

// valued that I need to characterize the local signals in the vincinity. 
float local_val = 0 ;
int   NvoxInVinc = (2*vinc+1)*(2*vinc+1)*(2*vinc+1);
double vec1[NvoxInVinc] ; 
for(int it = 0 ; it < NvoxInVinc ; it++) vec1[it] = 0 ; 
float  grad_stdev = 0; 
float  value_dist = 0; 


/// for estimation and out put of program process and how much longer it will take.
 int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead; 
 int running_index = 0 ; 
 int pref_ratio = 0 ;

if (sizeSlice * sizePhase * sizeRead > 32767) cout << " the number of voxels is bigger than the range of int the time estimation will be wrong " << endl; 

if ( do_masking == 1 ) {
	nvoxels_to_go_across = 0; 
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
	          if (*(nim_roi_data  +  nxy*islice + nx*ix  + iy  ) > 0 ) {
        		 nvoxels_to_go_across = nvoxels_to_go_across +1 ;
        	  }	 
           } 
	    }
	  }
}
 
 
 cout << " The number of voxels to go across = "<< nvoxels_to_go_across << endl ; 

/////////////////////////
////SMOOTHING LOOP  /////
/////////////////////////
 cout << " Big smoothing loop is beeing done now" << endl ; 



	for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
		 if ( !( !(*(nim_roi_data  +  nxy*iz + nx*ix  + iy  ) > 0) && (do_masking == 1)  ) ) {
         //if (iz==sizeSlice/2 && iy == sizePhase/2-4 && ix == sizeRead/2-4 ) { // debug loop open
         
           // this is to write out how many more voxels I have to go through. 
            running_index ++ ; 
            if ((running_index*100)/nvoxels_to_go_across != pref_ratio ) {
         	   cout << "\r "<<(running_index*100)/nvoxels_to_go_across <<  "% is done" << flush ; 
         	   pref_ratio = (running_index*100)/nvoxels_to_go_across ; 
            }
         
          // I am cooking in a clean kitchen. 
          // this cleaning might not be necesary, just to be on the save side
          *(gausweight_data  + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  *(smoothed_data    + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  NvoxInVinc = 0; 
		  local_val = *(nim_mask_data  + nxy*iz + nx*ix  + iy  ) ; 
		  
		  //examining the environment.
		  // and determining what the signal intensities are and what its distribution are 
		  for(int iz_i=max(0,iz-vinc); iz_i<=min(iz+vinc,sizeSlice-1); ++iz_i){
	        for(int iy_i=max(0,iy-vinc); iy_i<=min(iy+vinc,sizePhase-1); ++iy_i){
	      	 for(int ix_i=max(0,ix-vinc); ix_i<=min(ix+vinc,sizeRead-1);   ++ix_i){

		       vec1[NvoxInVinc] = (double) *(nim_mask_data  + nxy*iz_i + nx*ix_i  + iy_i  ) ; 
		       NvoxInVinc ++ ;

		  	 }	  
	        }
	       }
	       
		  // the standard deviation of the sinal valued in the vicinity, 
		  // this is necessary to normalice how many voxels are contributing to the local smoothing. 
		  //grad_stdev = (float )  gsl_stats_sd (vec1, 1, NvoxInVinc); 
		  grad_stdev = (float ) ren_stdev (vec1, NvoxInVinc); 
		
			for(int iz_i=max(0,iz-vinc); iz_i<=min(iz+vinc,sizeSlice-1); ++iz_i){
	    		for(int iy_i=max(0,iy-vinc); iy_i<=min(iy+vinc,sizePhase-1); ++iy_i){
	      			for(int ix_i=max(0,ix-vinc); ix_i<=min(ix+vinc,sizeRead-1); ++ix_i){
		  				dist_i     = dist((float)ix,(float)iy,(float)iz,(float)ix_i,(float)iy_i,(float)iz_i,dX,dY,dZ); 
		  				value_dist = fabs(  local_val - *(nim_mask_data  + nxy*iz_i + nx*ix_i  + iy_i)   );

		  					// *(debug_data       + nxy*iz_i + nx*ix_i  + iy_i  ) =  gaus(dist_i ,FWHM_val )/gaus(0,FWHM_val )  
		  					//													* gaus(value_dist,grad_stdev*0.1) /gaus(0,grad_stdev*0.1) ;  
		  					
		  				temp_wight_factor = gaus(dist_i ,FWHM_val ) *  gaus(value_dist,grad_stdev*selectivity)/gaus(0,grad_stdev*selectivity) ; 
		  				
		  				// The gaus data are important to avoid local scaling differences, when the kernel size changes. E.g. at edge of images.
		  				// this is a geometric parameter and only need to be calculated for one time point. 
		  				// this might be avoidable, if the gaus fucnction is better normaliced. 
		  			    *(gausweight_data  + nxy*iz + nx*ix  + iy  ) = 	 *(gausweight_data  + nxy*iz + nx*ix  + iy  ) 
		    														     + temp_wight_factor ;   			
			  		        
		  				for(int it=0; it<nrep; ++it){    // loop across lall time steps
		  					*(smoothed_data + nxyz *it  + nxy*iz + nx*ix  + iy  ) =   *(smoothed_data    + nxyz *it + nxy*iz   + nx*ix    + iy    ) 
		  					                                                        + *(nim_inputf_data  + nxyz *it + nxy*iz_i + nx*ix_i  + iy_i  ) * temp_wight_factor ;
		    			}	   			
			  		        
		            }	  
	      	    }
	       }
	     // scaling the signal intensity with the overall gaus leakage 
	     if (*(gausweight_data  + nxy*iz + nx*ix  + iy  ) > 0 ) {
	       	for(int it=0; it<nrep; ++it){  
	          *(smoothed_data   + nxyz *it   + nxy*iz + nx*ix  + iy  )  = *(smoothed_data  + nxyz *it   + nxy*iz + nx*ix  + iy  )/ *(gausweight_data  + nxy*iz + nx*ix  + iy  );
	        }
	     }
	     //if (*(nim_mask_data   +  nxy*iz + nx*ix  + iy  )  <= 0 )	     	*(smoothed_data    + nxy*iz + nx*ix  + iy  ) =  *(nim_inputf_data  + nxy*iz + nx*ix  + iy  ) ;
	    
	     //}//debug loop closed
	     }
        }
      }
    }

   cout << endl; 


// I am not sure if there is a case there masking makes sense? 
// I just leave it in. 

*/
//cout << " runing also until here  5.... " << endl; 
//cout << " slope " << smoothed->scl_slope << " " << nim_inputfi->scl_slope  << endl; 


  
    save_output_nifti(fout, "", nii_leak, true, use_outpath);
//  const char  *fout_2="debug.nii" ;
//  if( nifti_set_filenames(debug, fout_2 , 1, 1) ) return 1;
//  nifti_image_write( debug);

  return 0;
}


