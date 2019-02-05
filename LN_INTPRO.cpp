
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
      "LN_INTPRO: A program that does maximun and minimum intensity projections \n"
      "\n"
      "    This is very helpful for venous vessle imaging \n"
      "\n"
      "\n"
      "    basic usage: LN_INTPRO -image file.nii -min -direction 3 \n"
      "\n"
      "\n"
      "\n"
      "    options:\n"
      "\n"
      "       -help               : show this help\n"
      "       -image file.nii     : nii image that should be used for intensity projections \n"
      "       -max                : do maximum intentisy projection (do not combine with min) \n"
      "       -min                : do minimum intentisy projection (do not combine with max) \n"
      "       -direction 1        : direction in wich the dimention is collapst \n"
      "                               1 stand for x, 2 stand for y, and 3 stands for slice \n"
      "       -range              : optional parameter for \n"
      "                               the range of neigbouring voxles included\n"
      "                               default is all aslices \n"
      "\n"
      "\n"
      "    If there is a time series, the entire time domain is collapsed too.\n"
      "    E.g. the program chooses the minimum/maximmum value across time and slices\n"
      "\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   //nifti_image * nim_input=NULL;
   char        * fin_1=NULL;
   int          ac, is_max=0, is_min= 0;
   int 			is_direction=0, is_range=0;
   if( argc < 2 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-direction") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -direction\n");
            return 1;
         }
         is_direction = atoi(argv[ac]);  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-range") ) {
		  if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -range\n");
            return 1;
         }
         is_range = atoi(argv[ac]);  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-min") ) {
         is_min = 1;  // no string copy, just pointer assignment 
         cout << " I am doing MINIMUM intensity projections " << endl; 
      }
      else if( ! strcmp(argv[ac], "-max") ) {
         is_max = 1;  // no string copy, just pointer assignment 
         cout << " I am doing MAXIMUM intensity projections " << endl; 
      }
      else if( ! strcmp(argv[ac], "-image") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -timeseries\n");
            return 1;
         }
         fin_1 = argv[ac];  // no string copy, just pointer assignment 
      }
   }

   if( !fin_1  ) { fprintf(stderr, "** missing option '-timeseries'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_file_1i = nifti_image_read(fin_1, 1);
   if( !nim_file_1i ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin_1);
      return 2;
   }

   if( is_min+is_max !=1  ) { fprintf(stderr, "** you are supposed to do either maximum or minimum, ONE of each \n");  return 1; }
   
   
   if (is_direction-2 >1) {cout << " I don't understand this direction input " ; return 1 ; }
   
   cout << " the signal is collapsed along " ;
   if (is_direction==1) cout << " x-direction "<< endl ;
   if (is_direction==2) cout << " y-direction "<< endl ;
   if (is_direction==3) cout << " slice-direction "<< endl ;

   // get dimsions of input 
   int sizeSlice = nim_file_1i->nz ; 
   int sizePhase = nim_file_1i->nx ; 
   int sizeRead = nim_file_1i->ny ; 
   int nrep =  nim_file_1i->nt; 
   int nx =  nim_file_1i->nx;
   int nxy = nim_file_1i->nx * nim_file_1i->ny;
   int nxyz = nim_file_1i->nx * nim_file_1i->ny * nim_file_1i->nz;

   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 

   if (is_range > 0 ) cout << " I am ony looking for the Maximum/Minimum value in a proximity of " << is_range << " voxels " << endl; 
   if (is_range == 0 ) is_range = max (max(sizeSlice,sizePhase),sizeRead)  ;

   nifti_image * nim_file_1  	= nifti_copy_nim_info(nim_file_1i);
   nim_file_1->datatype = NIFTI_TYPE_FLOAT32;
   nim_file_1->nbyper = sizeof(float);
   nim_file_1->data = calloc(nim_file_1->nvox, nim_file_1->nbyper);
   float  *nim_file_1_data = (float *) nim_file_1->data;
   


  // if( !fout ) { fprintf(stderr, "-- no output requested \n"); return 0; }
   // assign nifti_image fname/iname pair, based on output filename
   //   (request to 'check' image and 'set_byte_order' here) 
  // if( nifti_set_filenames(nim_input, fout, 1, 1) ) return 1;

if ( nim_file_1i->datatype == NIFTI_TYPE_FLOAT32 ) {
  float  *nim_file_1i_data = (float *) nim_file_1i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int iz=0; iz<sizeSlice; ++iz){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_1_data  + nxyz *it +  nxy*iz + nx*ix  + iy  ) = (float) (*(nim_file_1i_data  + nxyz *it +  nxy*iz + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}  
  

if ( nim_file_1i->datatype == NIFTI_TYPE_INT16 ) {
  short  *nim_file_1i_data = (short *) nim_file_1i->data;
  	for(int it=0; it<nrep; ++it){  
	  for(int iz=0; iz<sizeSlice; ++iz){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
        		 *(nim_file_1_data  + nxyz *it +  nxy*iz + nx*ix  + iy  ) = (float) (*(nim_file_1i_data  + nxyz *it +  nxy*iz + nx*ix  + iy  )) ;	
           } 
	    }
	  }
	}
}    

    nifti_image * collaps_file  	= nifti_copy_nim_info(nim_file_1);
    collaps_file->nt 				= 1	; 
    collaps_file->nvox 				= nim_file_1->nvox / nrep  ; 
    collaps_file->datatype 			= NIFTI_TYPE_FLOAT32; 
    collaps_file->nbyper 			= sizeof(float);
    collaps_file->data 				= calloc(collaps_file->nvox, collaps_file->nbyper);
    float  *collaps_file_data 	    = (float *) collaps_file->data;

  
	double vec_file1[nrep]  ;

    /////////////////////////////////////////////////
    cout << " Starting with dimentionality collapse  = " << endl ; 
    /////////////////////////////////////////////////
    
    float extreme_val  = 0.0; 
    if (is_min==1) extreme_val = 100000000000.0;
    
if (is_direction == 3) {   
     for(int ix=0; ix<sizeRead; ++ix){
		 for(int iy=0; iy<sizePhase; ++iy){

           for(int iz=0; iz<sizeSlice; ++iz){ 	
			if (is_min==1) extreme_val = 100000000000.0;
		    if (is_max==1) extreme_val = 0.0;  	   
			for(int iz_i=max(iz-is_range,0); iz_i<min(iz+is_range,sizeSlice); ++iz_i){    
	         for(int it=0; it<nrep; ++it){ 
				
				if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz_i + nx*ix  + iy  ) < extreme_val  && *(nim_file_1_data  + nxyz*it + nxy*iz_i + nx*ix  + iy  ) > 0.0 ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz_i + nx*ix  + iy  ) ;
				} 
				if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz_i + nx*ix  + iy  ) > extreme_val  ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz_i + nx*ix  + iy  ) ;
				} 
								
			 } 
			}   
	        for(int it=0; it<nrep; ++it){ 
               if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) >= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ; 
               if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) <= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ;  
            }
	       }
	    }
	  }
}


if (is_direction == 1) {   
     for(int iz=0; iz<sizeSlice; ++iz){
		 for(int iy=0; iy<sizePhase; ++iy){

           for(int ix=0; ix<sizeRead; ++ix){ 	
			if (is_min==1) extreme_val = 100000000000.0;
		    if (is_max==1) extreme_val = 0.0;  	   
			for(int ix_i=max(ix-is_range,0); ix_i<min(ix+is_range,sizeRead); ++ix_i){    
	         for(int it=0; it<nrep; ++it){ 
				
				if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix_i  + iy  ) < extreme_val  && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix_i  + iy  ) > 0.0 ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix_i  + iy  ) ;
				} 
				if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix_i  + iy  ) > extreme_val  ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix_i  + iy  ) ;
				} 
								
			 } 
			}   
	        for(int it=0; it<nrep; ++it){ 
               if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) >= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ; 
               if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) <= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ;  
            }
	       }
	    }
	  }
}

if (is_direction == 2) {   
     for(int iz=0; iz<sizeSlice; ++iz){
		 for(int ix=0; ix<sizeRead; ++ix){

           for(int iy=0; iy<sizePhase; ++iy){ 	
			if (is_min==1) extreme_val = 100000000000.0;
		    if (is_max==1) extreme_val = 0.0;  	   
			for(int iy_i=max(iy-is_range,0); iy_i<min(iy+is_range,sizePhase); ++iy_i){    
	         for(int it=0; it<nrep; ++it){ 
				
				if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy_i  ) < extreme_val  && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy_i  ) > 0.0 ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy_i  ) ;
				} 
				if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy_i  ) > extreme_val  ) {
					extreme_val = *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy_i  ) ;
				} 
								
			 } 
			}   
	        for(int it=0; it<nrep; ++it){ 
               if( is_min==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) >= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ; 
               if( is_max==1 && *(nim_file_1_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) <= extreme_val ) *(collaps_file_data  + nxyz*it + nxy*iz + nx*ix  + iy  ) =  extreme_val ;  
            }
	       }
	    }
	  }
}

  string prefix_1 = "collaps_" ;
  string filename_1 = (string) (fin_1) ;
  string outfilename_1 = prefix_1+filename_1 ;
  cout << "writing colapsed file as = " << outfilename_1.c_str() << endl; // finfi is: char *
  const char  *fout_1=outfilename_1.c_str() ;
  if( nifti_set_filenames(collaps_file, fout_1 , 1, 1) ) return 1;
  nifti_image_write( collaps_file );

  
  

 // const char  *fout_6="kootrGM.nii" ;
 // if( nifti_set_filenames(GMkoord2, fout_6 , 1, 1) ) return 1;
 // nifti_image_write( GMkoord2 );

 // koord.autowrite("koordinaten.nii", wopts, &prot);
  return 0;
}



 
