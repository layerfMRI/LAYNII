
// Ausf¸hren mit . ./layers border_example_resized.nii brain_maskexample_resized.nii 0

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
//#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
using namespace std;

#define PI 3.14159265; 

//#include "utils.hpp"

int show_help( void )
{
   printf(
      "LN_CORREL2FILES: To estimate the voxel wise correlation of two timeseries\n"
      "\n"
      "    This program is motivated by Eli Merriam comparing in hunting down voxels that out of phase for VASO and BOLD      \n"
      "\n"
      "\n"
      "    basic usage: LN_CORREL2FILES -file1 file1.nii -file2 file2.nii \n"
      "\n"
      "\n"
      "\n"
      "    options:\n"
      "\n"
      "       -help               : show this help\n"
      "       -file1 			  : first time series\n"
      "       -file2  	          : second time series with should have the same dimensions as first time series\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   //nifti_image * nim_input=NULL;
   char        * fin_1=NULL, * fin_2=NULL ;
   int          ac, disp_float_eg=0;
   if( argc < 2 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-file1") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -file1\n");
            return 1;
         }
         fin_1 = argv[ac];  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-file2") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -file2\n");
            return 1;
         }
         fin_2 = argv[ac];  // no string copy, just pointer assignment 
      }
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !fin_1  ) { fprintf(stderr, "** missing option '-file1'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_file_1i = nifti_image_read(fin_1, 1);
   if( !nim_file_1i ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin_1);
      return 2;
   }
   
   if( !fin_2  ) { fprintf(stderr, "** missing option '-file2'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_file_2i = nifti_image_read(fin_2, 1);
   if( !nim_file_2i ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin_2);
      return 2;
   }
   

   // get dimsions of input 
   int sizeSlice = nim_file_1i->nz ; 
   int sizePhase = nim_file_1i->nx ; 
   int sizeRead = nim_file_1i->ny ; 
   int nrep =  nim_file_1i->nt; 
   int nx =  nim_file_1i->nx;
   int nxy = nim_file_1i->nx * nim_file_1i->ny;
   int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 


   nifti_image * nim_file_1  	= nifti_copy_nim_info(nim_file_1i);
   nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
   nim_file_1->nbyper = sizeof(float);
   nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
   float  *nim_file_1_data = (float *) nim_file_1->data;
   
   nifti_image * nim_file_2  	= nifti_copy_nim_info(nim_file_1i);
   nim_file_2->datatype = NIFTI_TYPE_FLOAT32;
   nim_file_2->nbyper = sizeof(float);
   nim_file_2->data = calloc(nim_file_2->nvox, nim_file_2->nbyper);
   float  *nim_file_2_data = (float *) nim_file_2->data;


  // if( !fout ) { fprintf(stderr, "-- no output requested \n"); return 0; }
   // assign nifti_image fname/iname pair, based on output filename
   //   (request to 'check' image and 'set_byte_order' here) 
  // if( nifti_set_filenames(nim_input, fout, 1, 1) ) return 1;

if ( nim_file_1i->datatype == NIFTI_TYPE_FLOAT32 ) {
  float  *nim_file_1i_data = (float *) nim_file_1i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_1_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_file_1i_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}  
  

if ( nim_file_1i->datatype == NIFTI_TYPE_INT16 ) {
  short  *nim_file_1i_data = (short *) nim_file_1i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_1_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_file_1i_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    

if ( nim_file_1i->datatype == NIFTI_TYPE_FLOAT32 ) {
  float  *nim_file_1i_data = (float *) nim_file_1i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_1_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_file_1i_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}  
  

if ( nim_file_2i->datatype == NIFTI_TYPE_INT16 ) {
  short  *nim_file_2i_data = (short *) nim_file_2i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_2_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_file_2i_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}   


if ( nim_file_2i->datatype == NIFTI_TYPE_FLOAT32 ) {
  float  *nim_file_2i_data = (float *) nim_file_2i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_2_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = (float) (*(nim_file_2i_data  + nxyz *it +  nxy*islice + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}  


    nifti_image * correl_file  		= nifti_copy_nim_info(nim_file_1);
    correl_file->nt 				= 1	; 
    correl_file->nvox 				= nim_file_1->nvox / nrep ; 
    correl_file->datatype 			= NIFTI_TYPE_FLOAT32; 
    correl_file->nbyper 			= sizeof(float);
    correl_file->data 				= calloc(correl_file->nvox, correl_file->nbyper);
    float  *correl_file_data 		= (float *) correl_file->data;

double vec_file1[nrep]  ;
double vec_file2[nrep]  ;

	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
	          	for(int it=0;  it<nrep; ++it){  
					vec_file1[it] = *(nim_file_1_data + nxyz*it + nxy*islice + nx*ix  + iy  );
	        		vec_file2[it] = *(nim_file_2_data + nxyz*it + nxy*islice + nx*ix  + iy  );
	            }
        		*(correl_file_data +  nxy*islice + nx*ix  + iy  ) =  gsl_stats_correlation(vec_file1,1,vec_file2,1,nrep)  ;	
           } 
	    }
	  }

 cout << " runing also until here  5.... " << endl; 

  string prefix = "correlated_" ;
  string filename_1 = (string) (fin_1) ;
  string outfilename = prefix+filename_1 ;
  
   cout << "writing as = " << outfilename.c_str() << endl; // finfi is: char *

  const char  *fout_1=outfilename.c_str() ;
  if( nifti_set_filenames(correl_file, fout_1 , 1, 1) ) return 1;
  nifti_image_write( correl_file );

 // const char  *fout_5="debug_ing.nii" ;
 // if( nifti_set_filenames(growfromWM0, fout_5 , 1, 1) ) return 1;
 // nifti_image_write( growfromWM0 );
  
 // const char  *fout_6="kootrGM.nii" ;
 // if( nifti_set_filenames(GMkoord2, fout_6 , 1, 1) ) return 1;
 // nifti_image_write( GMkoord2 );

 // koord.autowrite("koordinaten.nii", wopts, &prot);
  return 0;
}



 
