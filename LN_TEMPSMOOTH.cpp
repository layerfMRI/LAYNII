
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
//#include <gsl/gsl_statistics_double.h>
using namespace std;

#define PI 3.14159265; 


//#include "utils.hpp"

int show_help( void )
{
   printf(
      "LN_TEMPSMOOTH : Temporal filtering\n"
      "\n"
      "    This program is designet smooth data within the time domain ,\n"
      "    It is basically like a low pass filter (removing high frequencies),\n"
      "    I use a gouassian or box-car weigth function\n"
      "\n"
      "\n"
      "    basic usage: LN_TEMPSMOOTH -timeseries file.nii -gaus 1.0  \n"
      "    basic usage: LN_TEMPSMOOTH -timeseries file.nii -box 1  \n"
      "\n"
      "\n"
      "   This program now supports INT16, INT32 and FLOAT23 \n"

      "\n"
      "       -help                 : show this help\n"
      "       -timeseries           : nii file that with the series that should be smoothed \n"
      "                               only the first time point of this file is used.  \n"
      "       -gaus  value          : doing the smoothing with a box car, basically a trawiling window of averaging  \n"
      "                               specify the value of the gausian sice (float values) uin units of TR \n"
      "       -box   value          : doing the smoothing with a gausian  \n"
      "                               specify the value of the box sice (intager value) \n"
      "       -help                 : show this help\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   char       * fmaski=NULL, * fout=NULL, * finfi=NULL ;
   int          ac , do_gaus = 0 , do_box = 0, bFWHM_val=0   ;
   float 		gFWHM_val=0.0 ;
   if( argc < 3 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-gaus") ) {
        if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -gaus\n");
            return 1;
         }
         gFWHM_val = atof(argv[ac]);  // no string copy, just pointer assignment 
         do_gaus = 1; 
         cout << " I will do gaussian temporal smoothing " << endl;
      }
      else if( ! strcmp(argv[ac], "-box") ) {
        if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -box\n");
            return 1;
         }
         bFWHM_val = atoi(argv[ac]);  // no string copy, just pointer assignment 
         do_box = 1; 
          cout << " I will do box car like smoothing, this is like a running average. sliding window" << endl;
      }
      else if( ! strcmp(argv[ac], "-timeseries") ) {
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

   if( !finfi  ) { fprintf(stderr, "** missing option '-timeseries'\n");  return 1; }
   // read input dataset, including data 
   nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
   if( !nim_inputfi ) {
      fprintf(stderr,"** failed to read layer NIfTI image from '%s'\n", finfi);
      return 2;
   }
   

   
   if (do_box + do_gaus !=1) {
   	  cout << " I don't know which kind of smoothing I should use (gaussian or box car), please decide which one it should be " << endl;
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
   float dX =  nim_inputfi->pixdim[1] ; 
   float dY =  nim_inputfi->pixdim[2] ; 
   float dZ =  nim_inputfi->pixdim[3] ; 
   
   
// if you are running the smoothing in 2D, it will still go thought he entire pipeline. 
// the only difference is that the weights in a certain direction are suppressed
// doing it in 2Dim, will not speed up the program
   
   nifti_image * nim_inputf  	= nifti_copy_nim_info(nim_inputfi);
   nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
   nim_inputf->nbyper = sizeof(float);
   nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
   float  *nim_inputf_data = (float *) nim_inputf->data;
   
   
   
   /////////////////////////////////////////////////////////////////////////
   /////////  fixing potential problems with different input datatypes /////
   /////////  here, I am loading them in their native datatype /////////////
   /////////  and translate them to the datatime I like best  //////////////
   /////////////////////////////////////////////////////////////////////////

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


	// write out some stuff that might be good to know, if you want to debug
   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 
   cout << " Voxel size    " <<  dX << " x " <<  dY << " x "  <<  dZ  << endl; 

	cout << " datatye 1 = " << nim_inputf->datatype << endl;

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



//float dist (float x1, float y1, float z1, float x2, float y2, float z2,float dX, float dY, float dZ) ; 
float gaus (float distance, float sigma) ;

cout << "debug  2 " << endl; 


/////////////////////////////////////////////////////////////////////////////
////SMOOTHING LOOP  for the case you chose gaussian smoothing /////
/////////////////////////////////////////////////////////////////////////////
if (do_gaus) { /// 
	
	
//float kernal_size = 10; // corresponds to one voxel sice. 
int vinc = max(1.,2. * gFWHM_val/dX ); // if voxel is too far away, I ignore it. 
float dist_i = 0.;
cout << " vinc " <<  vinc<<  endl; 
cout << " FWHM_val " <<  gFWHM_val<<  endl; 


/// for estimation and out put of program process and how much longer it will take.
int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead; 
int running_index = 0 ; 
int pref_ratio = 0 ;


///////////////////////////////////
////SMOOTHING LOOP  /////
///////////////////////////////////
//cout << " DEBUG " <<   dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl; 

 cout << " smoothing with gaus "  << endl ; 

	for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
			
			
			 running_index ++ ; 
            if ((running_index*100)/nvoxels_to_go_across != pref_ratio ) {
         	   cout << "\r "<<(running_index*100)/nvoxels_to_go_across <<  "% is done " << flush ; 
         	   pref_ratio = (running_index*100)/nvoxels_to_go_across ; 
            }
			

			
			
		 for(int it=0; it<nrep; ++it){	
           *(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  //*(smoothed_data    + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  
	      			for(int it_i=max(0,it-vinc); it_i<min(it+vinc+1,nrep); ++it_i){
	      			  if ( *(nim_inputf_data + nxyz *it_i + nxy*iz + nx*ix  + iy  ) != 0 ) {
		  				dist_i = abs (it-it_i ); 
		  				//cout << "debug  4 " <<  gaus(dist_i ,FWHM_val ) <<   endl; 
		  			    //cout << "debug  5 " <<  dist_i  <<   endl; 
						//if ( *(nim_input_data   +  nxy*iz + nx*ix  + iy  )  == 3 ) cout << "debug  4b " << endl; 

							//dummy = *(layer_data  + nxy*iz_i + nx*ix_i  + iy_i  ); 
		  					*(smoothed_data   + nxyz *it + nxy*iz + nx*ix  + iy  ) = *(smoothed_data   + nxyz *it + nxy*iz + nx*ix  + iy  ) + *(nim_inputf_data + nxyz *it_i + nxy*iz + nx*ix  + iy  ) * gaus(dist_i ,gFWHM_val ) ;
		    				*(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  ) = *(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  ) + gaus(dist_i ,gFWHM_val ) ; 
							
			  		  }	
			  	  }	

	     
        }
      }
    }
  }  

 cout  << endl ; 
  
}

/////////////////////////////////////////////////////////////////////////////
//// CLOSED SMOOTHING LOOP  for the case you chose gaussian smoothing CLOSED /////
/////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////
////SMOOTHING LOOP  for the case you chose box car smoothing /////
/////////////////////////////////////////////////////////////////////////////
if (do_box) { /// 
	
	
//float kernal_size = 10; // corresponds to one voxel sice. 
int vinc = bFWHM_val; // if voxel is too far away, I ignore it. 
cout << " vinc " <<  vinc<<  endl; 
float dist_i = 0.;

/// for estimation and out put of program process and how much longer it will take.
int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead; 
int running_index = 0 ; 
int pref_ratio = 0 ;


///////////////////////////////////
////SMOOTHING LOOP  /////
///////////////////////////////////
//cout << " DEBUG " <<   dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl; 

 cout << " smoothing with box car filter "  << endl ; 

	for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
			
			
			 running_index ++ ; 
            if ((running_index*100)/nvoxels_to_go_across != pref_ratio ) {
         	   cout << "\r "<<(running_index*100)/nvoxels_to_go_across <<  "% is done " << flush ; 
         	   pref_ratio = (running_index*100)/nvoxels_to_go_across ; 
            }
			

			
			
		 for(int it=0; it<nrep; ++it){	
           *(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  //*(smoothed_data    + nxy*iz + nx*ix  + iy  )  = 0 ; 
		  
	      			for(int it_i=max(0,it-vinc); it_i<min(it+vinc+1,nrep); ++it_i){
	      			  if ( *(nim_inputf_data + nxyz *it_i + nxy*iz + nx*ix  + iy  ) != 0 ) {
		  				// No distance here.... just averaging
		  				//dist_i = abs (it-it_i ); 

		  					*(smoothed_data   + nxyz *it + nxy*iz + nx*ix  + iy  ) = *(smoothed_data   + nxyz *it + nxy*iz + nx*ix  + iy  ) + *(nim_inputf_data + nxyz *it_i + nxy*iz + nx*ix  + iy  ) ;
		    				*(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  ) = *(gausweight_data + nxyz *it + nxy*iz + nx*ix  + iy  ) + 1 ; 
							
			  		  }	
			  	  }	

	     
        }
      }
    }
  }  

 cout  << endl ; 
  
}

/////////////////////////////////////////////////////////////////////////////
//// CLOSED SMOOTHING LOOP  for the case you chose box smoothing CLOSED /////
/////////////////////////////////////////////////////////////////////////////




 ///////////////////////////////////
 //// correcting for edge error  /////
 ///////////////////////////////////
 for(int it=0; it<nrep; ++it){  
  for(int iz=0; iz<sizeSlice; ++iz){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
		  if ( *(gausweight_data  + nxyz *it + nxy*iz + nx*ix  + iy  ) != 0 ) {
		  	*(smoothed_data    + nxyz *it + nxy*iz + nx*ix  + iy  ) = *(smoothed_data   + nxyz *it  + nxy*iz + nx*ix  + iy  ) / *(gausweight_data  + nxyz *it + nxy*iz + nx*ix  + iy  );							
          }
        }
      }
    }
 }
  

//cout << " runing also until here  5.... " << endl; 
//cout << " slope " << smoothed->scl_slope << " " << nim_inputfi->scl_slope  << endl; 

smoothed->scl_slope =  nim_inputfi->scl_slope  ;


if (nim_inputfi->scl_inter != 0 ){
cout << " ############################################################# " << endl; 
cout << " #############   WARNING   WANRING   WANRING  ################ " << endl; 
cout << " ########   the NIFTI scale factor is asymmetric  ############ " << endl; 
cout << " ########   Why would you do such a thing????     ############ " << endl; 
cout << " #############   WARNING   WANRING   WANRING  ################ " << endl; 
cout << " ############################################################# " << endl; 
}

  string prefix = "smoothed_" ;
  string filename = (string) (finfi) ;
  string outfilename = prefix+filename ;
  
  cout << "writing as = " << outfilename.c_str() << endl; // finfi is: char *

  const char  *fout_1=outfilename.c_str() ;
  if( nifti_set_filenames(smoothed, fout_1 , 1, 1) ) return 1;
  nifti_image_write( smoothed );
  
//  const char  *fout_2="debug.nii" ;
//  if( nifti_set_filenames(debug, fout_2 , 1, 1) ) return 1;
//  nifti_image_write( debug);

  return 0;
}


//  float dist (float x1, float y1, float z1, float x2, float y2, float z2, float dX, float dY, float dZ ) {
//    return sqrt((x1-x2)*(x1-x2)*dX*dX+(y1-y2)*(y1-y2)*dY*dY+(z1-z2)*(z1-z2)*dZ*dZ);
//  }

  float gaus (float distance, float sigma) {
    return 1./(sigma*sqrt(2.*3.141592))*exp (-0.5*distance*distance/(sigma*sigma));
  }

