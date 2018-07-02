
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
      "LN_RAGRUG: generating a file where every voxel has a specific value that is indicative of its position\n"
      "\n"
      "    This program is inspired by Kendrick Kay's way of visualizing which position at the surface comes from which voxel ,\n"
      "    see Fig. 4 of https://www.biorxiv.org/content/early/2018/06/03/337667 ,\n"
      "\n"
      "    basic usage: LN_RAGRUG -input any_file.nii  \n"
      "\n"
      "\n"
      "   This program now supports INT16, INT32 and FLOAT23 \n"

      "\n"
      "       -help               	: show this help\n"
      "       -input 		     	: nii file. this program will use the dimension of this file to generate a Rag Rug file accordingly /n"


      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   char       * fout=NULL, * finfi=NULL ;
   int          ac, direction_i=0; 
   float 		FWHM_val=0 ;
   if( argc < 3 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-input") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         finfi = argv[ac];  // no string copy, just pointer assignment 
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
   
   
      // get dimsions of input 
   int sizeSlice = nim_inputfi->nz ; 
   int sizePhase = nim_inputfi->nx ; 
   int sizeRead = nim_inputfi->ny ; 
   int nrep =  nim_inputfi->nt; 
   int nx =  nim_inputfi->nx;
   int nxy = nim_inputfi->nx * nim_inputfi->ny;
   int nxyz = nim_inputfi->nx * nim_inputfi->ny * nim_inputfi->nz;
   float dX =  1;//nim_inputfi->pixdim[1] ; 
   float dY =  1;//nim_inputfi->pixdim[2] ; 
   float dZ =  1;//nim_inputfi->pixdim[3] ; 
   
   

      
   //nim_mask->datatype = NIFTI_TYPE_FLOAT32;
   //nim_mask->nbyper = sizeof(float);
   //nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
   
   
   nifti_image * nim_inputf  	= nifti_copy_nim_info(nim_inputfi);
   nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
   nim_inputf->nbyper = sizeof(float);
   nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
   float  *nim_inputf_data = (float *) nim_inputf->data;
   

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


   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 
   cout << " Voxel size    " <<  dX << " x " <<  dY << " x "  <<  dZ  << endl; 

	
	cout << " datatye 1 = " << nim_inputf->datatype << endl;

///////////////////////////////////
////   MAKE allocating necessary files  /////
///////////////////////////////////

    
    nifti_image * ragrug  	= nifti_copy_nim_info(nim_inputf);
    nifti_image * coord  	= nifti_copy_nim_info(nim_inputf);

    ragrug->nt 	= 1	; 
	coord->nt   = 3 ; // for three spatial dimensions
	//nifti_update_dims_from_array(coord) ; 
	coord->nvox = nim_inputfi->nvox * 3 ; 

    ragrug->datatype 		= NIFTI_TYPE_INT16; 
	coord->datatype 		= NIFTI_TYPE_INT16;

    ragrug->nbyper 			= sizeof(short);
    coord->nbyper 		    = sizeof(short);

cout << " ragrug->nvox " << ragrug->nvox << endl; 
cout << " coord ->nvox " << coord->nvox << endl; 
cout << " nim_inputfi  " << nim_inputfi->nvox << endl; 
cout << " nxyz         " << nxyz << endl; 


    ragrug->data = calloc(ragrug->nvox, ragrug->nbyper);
    coord->data  = calloc( coord->nvox,  coord->nbyper);

    short  *ragrug_data = (short *) ragrug->data;
    short  *coord_data  = (short *) coord->data;    


//cout << "debug  2 " << endl; 


///////////////////////////////////
////SMOOTHING LOOP  /////
///////////////////////////////////
//cout << " DEBUG " <<   dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl; 

 cout << " filling nii with spatial values " << endl ; 

  for(int iz=0; iz<sizeSlice; ++iz){  
     for(int iy=0; iy<sizePhase; ++iy){
       for(int ix=0; ix<sizeRead; ++ix){
	  					*(ragrug_data  + nxy*iz + nx*ix  + iy  ) = 0;							
       }
     }
   }

	for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
		  //*(smoothed_data    + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  		
		    	*(coord_data + nxyz *0 + nxy*iz + nx*ix  + iy  ) = ix ; 
		    	*(coord_data + nxyz *1 + nxy*iz + nx*ix  + iy  ) = iy ; 
		    	*(coord_data + nxyz *2 + nxy*iz + nx*ix  + iy  ) = iz ; 
		    	
		    	if (ix%2 ==0 ) *(ragrug_data + nxy*iz + nx*ix  + iy  ) = *(ragrug_data + nxy*iz + nx*ix  + iy  ) + 4  ; 
		    	if (iy%2 ==0 ) *(ragrug_data + nxy*iz + nx*ix  + iy  ) = *(ragrug_data + nxy*iz + nx*ix  + iy  ) + 2  ; 
		    	if (iz%2 ==0 ) *(ragrug_data + nxy*iz + nx*ix  + iy  ) = *(ragrug_data + nxy*iz + nx*ix  + iy  ) + 1  ; 
		    	
        }
      }
    }




///////////////////////////////////
//// correcting for edge error  /////
///////////////////////////////////

//   for(int iz=0; iz<sizeSlice; ++iz){  
//      for(int iy=0; iy<sizePhase; ++iy){
//        for(int ix=0; ix<sizeRead-0; ++ix){
//		  					*(smoothed_data    + nxy*iz + nx*ix  + iy  ) = *(smoothed_data    + nxy*iz + nx*ix  + iy  ) / *(gausweight_data  + nxy*iz + nx*ix  + iy  );							
//        }
//      }
//    }

// cout << " runing also until here  5.... " << endl; 


ragrug->scl_slope = 1 ;
coord->scl_slope = 1 ;


  const char  *fout_1="ragrug.nii" ;
  if( nifti_set_filenames(ragrug, fout_1 , 1, 1) ) return 1;
  nifti_image_write( ragrug );

  const char  *fout_2="coordinates.nii" ;
  if( nifti_set_filenames(coord, fout_2 , 1, 1) ) return 1;
  nifti_image_write( coord );



  return 0;
}





