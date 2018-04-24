
#include <stdio.h>
#include "nifti2_io.h"
//#include "nifti2.h"
//#include "nifti1.h"
//#include "nifticdf.h"
//#include "nifti_tool.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
using namespace std;

#define PI 3.14159265; 


//#include "utils.hpp"

int show_help( void )
{
   printf(
      "LN_LAYER_SMOOTH : Layering algorithm based on iterative smothing\n"
      "\n"
      "    This program is designet smooth data within layer or columns ,\n"
      "    In order to avoid smoothing across masks a crawler smoothed only across connected voxels ,\n"
      "\n"
      "    basic usage: LN_LAYER_SMOOTH -layer_file equi_dist_layers.nii -input activity_map.nii -FWHM 1 \n"
      "\n"
      "\n"
      "   This program now supports INT16, INT32 and FLOAT23 \n"

      "\n"
      "       -help               	: show this help\n"
      "       -layer_file     		: nii file that contains layer or column masks \n"
      "       -input 		     	: nii file that should be smoothed. it should have same dimentions as layer file\n"
      "       -FWHM 		     	: the amount of smoothing in mm\n"
      "       -twodim 		     	: optional argument to do smoothing in 2 Dim only \n"
      "       -mask 		     	: optional argument to mask activity outside of layers \n"


      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   char       * fmaski=NULL, * fout=NULL, * finfi=NULL ;
   int          ac, twodim=0, do_masking=0; 
   float 		FWHM_val=0 ;
   if( argc < 3 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-layer_file") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -layer_file\n");
            return 1;
         }
         fmaski = argv[ac];  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-FWHM") ) {
        if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -FWHM\n");
            return 1;
         }
         FWHM_val = atof(argv[ac]);  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-input") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         finfi = argv[ac];  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-twodim") ) {
         twodim = 1;
         cout << " I will do smoothing only in 2D"  << endl; 
      }
     else if( ! strcmp(argv[ac], "-mask") ) {
         do_masking = 1;
         cout << " I will set every to zero outside the layers"  << endl; 
      }
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !finfi  ) { fprintf(stderr, "** missing option '-input'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
   if( !nim_inputfi ) {
      fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", finfi);
      return 2;
   }
   
   if( !fmaski  ) { fprintf(stderr, "** missing option '-layer_file'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_maski = nifti_image_read(fmaski, 1);
   if( !nim_maski ) {
      fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", fmaski);
      return 2;
   }
   
      // get dimsions of input 
   int sizeSlice = nim_maski->nz ; 
   int sizePhase = nim_maski->nx ; 
   int sizeRead = nim_maski->ny ; 
   int nrep =  nim_maski->nt; 
   int nx =  nim_maski->nx;
   int nxy = nim_maski->nx * nim_maski->ny;
   int nxyz = nim_maski->nx * nim_maski->ny * nim_maski->nz;
   float dX =  nim_maski->pixdim[1] ; 
   float dY =  nim_maski->pixdim[2] ; 
   float dZ =  nim_maski->pixdim[3] ; 
   
   
   if  (twodim == 1) dZ = 1000 * dZ ; 
   
      
   //nim_mask->datatype = NIFTI_TYPE_FLOAT32;
   //nim_mask->nbyper = sizeof(float);
   //nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
   
   
   nifti_image * nim_inputf  	= nifti_copy_nim_info(nim_inputfi);
   nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
   nim_inputf->nbyper = sizeof(float);
   nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
   float  *nim_inputf_data = (float *) nim_inputf->data;
   
   nifti_image * nim_mask  	= nifti_copy_nim_info(nim_inputfi);
   nim_mask->datatype = NIFTI_TYPE_INT32;
   nim_mask->nbyper = sizeof(int);
   nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
   int  *nim_mask_data = (int *) nim_mask->data;
   /////////////////////////////////////////////
   /////////  fixing potential problems with different input datatypes  //////////////
   /////////////////////////////////////////////

if ( nim_inputfi->datatype == NIFTI_TYPE_FLOAT32 ) {
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
  

if ( nim_inputfi->datatype == NIFTI_TYPE_INT16 ) {
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


if ( nim_inputfi->datatype == NIFTI_TYPE_INT32 ) {
  int  *nim_inputfi_data = (int *) nim_inputfi->data;
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

if ( nim_maski->datatype == NIFTI_TYPE_FLOAT32 ) {
  float  *nim_maski_data = (float *) nim_maski->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (int) (*(nim_maski_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    
  
if ( nim_maski->datatype == NIFTI_TYPE_INT16 ) {
  short  *nim_maski_data = (short *) nim_maski->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (int) (*(nim_maski_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    

if ( nim_maski->datatype == NIFTI_TYPE_INT32 ) {
  int  *nim_maski_data = (int *) nim_maski->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_mask_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (int) (*(nim_maski_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    
	

   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 
   cout << " Voxel size    " <<  dX << " x " <<  dY << " x "  <<  dZ  << endl; 

	
	cout << " datatye 1 = " << nim_inputf->datatype << endl;
    cout << " datatye 2 = " << nim_mask ->datatype << endl;

///////////////////////////////////
////   MAKE allocating necessary files  /////
///////////////////////////////////

    
    nifti_image * smoothed  	= nifti_copy_nim_info(nim_inputf);
    nifti_image * gausweight  	= nifti_copy_nim_info(nim_inputf);
//    nifti_image * layer  		= nifti_copy_nim_info(nim_input);
//    nifti_image * leak_layer  	= nifti_copy_nim_info(nim_input);


    smoothed->datatype 		= NIFTI_TYPE_FLOAT32; 
	gausweight->datatype 	= NIFTI_TYPE_FLOAT32;
//	layer->datatype 		= NIFTI_TYPE_FLOAT32;
//	leak_layer->datatype 	= NIFTI_TYPE_FLOAT32;

    smoothed->nbyper 		= sizeof(float);
	gausweight->nbyper 		= sizeof(float);
//	layer->nbyper 			= sizeof(float);
//	leak_layer->nbyper 		= sizeof(float);

    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
//    layer->data = calloc(layer->nvox, layer->nbyper);
//    leak_layer->data = calloc(leak_layer->nvox, leak_layer->nbyper);

    float  *smoothed_data = (float *) smoothed->data;
    float  *gausweight_data = (float *) gausweight->data;    
//    float  *layer_data = (float *) layer->data;
//    float  *leak_layer_data = (float *) leak_layer->data;



float dist (float x1, float y1, float z1, float x2, float y2, float z2,float dX, float dY, float dZ) ; 
float gaus (float distance, float sigma) ;

cout << "debug  2 " << endl; 



//float kernal_size = 10; // corresponds to one voxel sice. 
int vinc = max(1.,2. * FWHM_val/dX ); // if voxel is too far away, I ignore it. 
float dist_i = 0.;
cout << " vinc " <<  vinc<<  endl; 
cout << " FWHM_val " <<  FWHM_val<<  endl; 

///////////////////////////////////
////Finding number of layers  /////
///////////////////////////////////
 int layernumber = 0 ;

	for(int iz=0; iz<sizeSlice; ++iz){
  		for(int iy=0; iy<sizePhase; ++iy){
       		for(int ix=0; ix<sizeRead; ++ix){
	 		 if (*(nim_mask_data   +  nxy*iz + nx*ix  + iy  ) > layernumber)  layernumber = *(nim_mask_data   +  nxy*iz + nx*ix  + iy) ; 
	  		  
      		} 
   		} 
   	}
   	
cout << " There are  " <<  layernumber<< " layers/masks to smooth within  " << endl; 


///////////////////////////////////
////making it intager   /////
///////////////////////////////////

///////////////////////////////////
////SMOOTHING LOOP  /////
///////////////////////////////////
//cout << " DEBUG " <<   dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl; 

 cout << " smoothing in layer   " << flush ; 

 for(int layernumber_i=1; layernumber_i<=layernumber; ++layernumber_i){  
 cout << "  " <<  layernumber_i << flush ; 

	for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
          *(gausweight_data  + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  //*(smoothed_data    + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  
	     if (*(nim_mask_data   +  nxy*iz + nx*ix  + iy  )  == layernumber_i ){
		
			for(int iz_i=max(0,iz-vinc); iz_i<min(iz+vinc+1,sizeRead); ++iz_i){
	    		for(int iy_i=max(0,iy-vinc); iy_i<min(iy+vinc+1,sizePhase); ++iy_i){
	      			for(int ix_i=max(0,ix-vinc); ix_i<min(ix+vinc+1,sizeRead); ++ix_i){
	      			  if (*(nim_mask_data   +  nxy*iz_i + nx*ix_i  + iy_i  )  == layernumber_i ){
		  				dist_i = dist((float)ix,(float)iy,(float)iz,(float)ix_i,(float)iy_i,(float)iz_i,dX,dY,dZ); 
		  				//cout << "debug  4 " <<  gaus(dist_i ,FWHM_val ) <<   endl; 
		  			    //cout << "debug  5 " <<  dist_i  <<   endl; 

						//if ( *(nim_input_data   +  nxy*iz + nx*ix  + iy  )  == 3 ) cout << "debug  4b " << endl; 

							//dummy = *(layer_data  + nxy*iz_i + nx*ix_i  + iy_i  ); 
		  					*(smoothed_data    + nxy*iz + nx*ix  + iy  ) = *(smoothed_data    + nxy*iz + nx*ix  + iy  ) + *(nim_inputf_data  + nxy*iz_i + nx*ix_i  + iy_i  ) * gaus(dist_i ,FWHM_val ) ;
		    				*(gausweight_data  + nxy*iz + nx*ix  + iy  ) = *(gausweight_data  + nxy*iz + nx*ix  + iy  ) + gaus(dist_i ,FWHM_val ) ; 
							
			  			
			  		  }	
		            }	  
	      	    }
	       }

	       if (*(gausweight_data  + nxy*iz + nx*ix  + iy  ) > 0 ) *(smoothed_data    + nxy*iz + nx*ix  + iy  )  = *(smoothed_data    + nxy*iz + nx*ix  + iy  )/ *(gausweight_data  + nxy*iz + nx*ix  + iy  );
	     }
	     
	     if (*(nim_mask_data   +  nxy*iz + nx*ix  + iy  )  <= 0 )	     	*(smoothed_data    + nxy*iz + nx*ix  + iy  ) =  *(nim_inputf_data  + nxy*iz + nx*ix  + iy  ) ;
	     
	     
	     
        }
      }
    }
}

 cout << endl; 


///////////////////////////////////
//// masking if it is it wanted   /////
///////////////////////////////////

if ( do_masking == 1 ) {
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
	          if (*(nim_mask_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) == 0 ) {
        		 *(smoothed_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = 0 ;	
        	  }	 
           } 
	    }
	  }
	}
} 

 cout << " runing also until here  5.... " << endl; 


     // output file name       
//  const char  *fout_4="leaky_layers.nii" ;
//  if( nifti_set_filenames(leak_layer, fout_4 , 1, 1) ) return 1;
//  nifti_image_write( leak_layer );

//  const char  *fout_5="input_file.nii" ;
//  if( nifti_set_filenames(nim_inputf, fout_5 , 1, 1) ) return 1;
//  nifti_image_write( nim_inputf );
  
///  const char  *fout_2="mask.nii" ;
//  if( nifti_set_filenames(nim_mask, fout_2 , 1, 1) ) return 1;
//  nifti_image_write( nim_mask );
  
  string prefix = "smoothed_" ;
  string filename = (string) (finfi) ;
  string outfilename = prefix+filename ;
  
   cout << "writing as = " << outfilename.c_str() << endl; // finfi is: char *

  const char  *fout_1=outfilename.c_str() ;
  if( nifti_set_filenames(smoothed, fout_1 , 1, 1) ) return 1;
  nifti_image_write( smoothed );
  
//  const char  *fout_1="layer.nii" ;
//  if( nifti_set_filenames(layer, fout_1 , 1, 1) ) return 1;
//  nifti_image_write( layer );





  return 0;
}




  float dist (float x1, float y1, float z1, float x2, float y2, float z2, float dX, float dY, float dZ ) {
    return sqrt((x1-x2)*(x1-x2)*dX*dX+(y1-y2)*(y1-y2)*dY*dY+(z1-z2)*(z1-z2)*dZ*dZ);
  }


  float gaus (float distance, float sigma) {
    return 1./(sigma*sqrt(2.*3.141592))*exp (-0.5*distance*distance/(sigma*sigma));
  }

