
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
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
using namespace std;

#define PI 3.14159265; 


#define PI 3.14159265;

//#include "utils.hpp"

int show_help( void )
{
   printf(
      "LN_NOISEME: short exmample of LAyering\n"
      "\n"
      "    This program is to demonstrate how to read a NIfTI-2 dataset,\n"
      "    set output filenames and write a NIfTI-2 dataset, all via the\n"
      "    standard NIfTI C library.\n"
      "\n"
      "    basic usage: LN_GROW_LAYERS -rim rim.nii \n"
      "\n"
      "\n"
      "   note that the rim.nii file always needs to be in datatype SHORT  \n"
      "\n"
      "    options:\n"
      "\n"
      "       -help               : show this help\n"
      "       -disp_float_example : show some voxel's data\n"
      "       -rim  border      : specify input dataset\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   nifti_image * nim_input=NULL;
   char        * fin=NULL, * fout=NULL;
   int          ac, disp_float_eg=0;
   if( argc < 2 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-disp_float_example") ) {
         disp_float_eg = 1;
      }
      else if( ! strcmp(argv[ac], "-rim") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         fin = argv[ac];  // no string copy, just pointer assignment 
      }
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !fin  ) { fprintf(stderr, "** missing option '-rim'\n");  return 1; }
   // read input dataset, including data 
   nim_input = nifti_image_read(fin, 1);
   if( !nim_input ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin);
      return 2;
   }
   

   // get dimsions of input 
   int sizeSlice = nim_input->nz ; 
   int sizePhase = nim_input->nx ; 
   int sizeRead = nim_input->ny ; 
   int nrep =  nim_input->nt; 
   int nx =  nim_input->nx;
   int nxy = nim_input->nx * nim_input->ny;
   int nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 


  // if( !fout ) { fprintf(stderr, "-- no output requested \n"); return 0; }
   // assign nifti_image fname/iname pair, based on output filename
   //   (request to 'check' image and 'set_byte_order' here) 
  // if( nifti_set_filenames(nim_input, fout, 1, 1) ) return 1;

	if ( nim_input->datatype != 4 ) {
		 //nim_input->datatype = NIFTI_TYPE_INT16 ;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
		cout << " YOU GAVE ME THE WRONG DATATYPE" << endl;
	}


	// get access to data of nim_input 
    short  *nim_input_data = (short *) nim_input->data;
    
    

    //nifti_brick_list   NB_orig, NB_select;
    //nifti_image      * nim_orig, * nim_select;
    //const int64_t                blist[5] = { 0, 0, 0, 0, 0 };
 
    //nim_orig   = nifti_image_read_bricks("rim.nii", 0, NULL,  &NB_orig);
    //nim_select = nifti_image_read_bricks("rim.nii", 5, blist, &NB_select);
    //update_nifti_image_for_brick_list( nim_orig,  1 );
    
   // allocating an additional nii 
   //const int64_t blist_2[1] = { 0};
   //nifti_brick_list NB_select_2;
   //nifti_image * growfromWM = nifti_image_read_bricks("rim.nii", 0, NULL, &NB_select_2) ;
   
    nifti_image * growfromWM0  = nifti_copy_nim_info(nim_input);
	growfromWM0->datatype = NIFTI_TYPE_FLOAT32;
	growfromWM0->nbyper = sizeof(float);
    growfromWM0->data = calloc(growfromWM0->nvox, growfromWM0->nbyper);
    float  *growfromWM0_data = (float *) growfromWM0->data;
   
   	nifti_image * growfromWM1  = nifti_copy_nim_info(nim_input);
	growfromWM1->datatype = NIFTI_TYPE_FLOAT32;
	growfromWM1->nbyper = sizeof(float);
    growfromWM1->data = calloc(growfromWM1->nvox, growfromWM1->nbyper);
    float  *growfromWM1_data = (float *) growfromWM1->data;
   
   //nifti_image * growfromWM1 = nifti_image_read(fin, 1);
   //float  *growfromWM1_data = (float *) growfromWM1->data;
   //nifti_image * growfromWM0 = nifti_image_read(fin, 1);
   //float  *growfromWM0_data = (float *) growfromWM0->data;
   //growfromWM->dim[4] = 1 ;  
  // nifti_update_dims_from_array(growfromWM);   // changing according sizes nt etc. 

    
    nifti_image * WMkoord0  = nifti_copy_nim_info(nim_input);
    nifti_image * WMkoord1  = nifti_copy_nim_info(nim_input);
    nifti_image * WMkoord2  = nifti_copy_nim_info(nim_input);
    nifti_image * WMkoord3  = nifti_copy_nim_info(nim_input);
    
    WMkoord0->datatype = NIFTI_TYPE_INT32; 
	WMkoord1->datatype = NIFTI_TYPE_INT32;
	WMkoord2->datatype = NIFTI_TYPE_INT32;
	WMkoord3->datatype = NIFTI_TYPE_INT32;

    WMkoord0->nbyper = sizeof(int);
	WMkoord1->nbyper = sizeof(int);
	WMkoord2->nbyper = sizeof(int);
	WMkoord3->nbyper = sizeof(int);

    WMkoord0->data = calloc(WMkoord0->nvox, WMkoord0->nbyper);
    WMkoord1->data = calloc(WMkoord1->nvox, WMkoord1->nbyper);
    WMkoord2->data = calloc(WMkoord2->nvox, WMkoord2->nbyper);
    WMkoord3->data = calloc(WMkoord3->nvox, WMkoord3->nbyper);

    int  *WMkoord0_data = (int *) WMkoord0->data;
    int  *WMkoord1_data = (int *) WMkoord1->data;
    int  *WMkoord2_data = (int *) WMkoord2->data;
    int  *WMkoord3_data = (int *) WMkoord3->data;


    nifti_image * growfromGM0  = nifti_copy_nim_info(nim_input);
	growfromGM0->datatype = NIFTI_TYPE_FLOAT32;
	growfromGM0->nbyper = sizeof(float);
    growfromGM0->data = calloc(growfromGM0->nvox, growfromGM0->nbyper);
    float  *growfromGM0_data = (float *) growfromGM0->data;
   
    nifti_image * growfromGM1  = nifti_copy_nim_info(nim_input);
	growfromGM1->datatype = NIFTI_TYPE_FLOAT32;
	growfromGM1->nbyper = sizeof(float);
    growfromGM1->data = calloc(growfromGM1->nvox, growfromGM1->nbyper);
    float  *growfromGM1_data = (float *) growfromGM1->data;

   //nifti_image * growfromGM1 = nifti_image_read(fin, 1);
   //float  *growfromGM1_data = (float *) growfromGM1->data;
   //nifti_image * growfromGM0 = nifti_image_read(fin, 1);
   //float  *growfromGM0_data = (float *) growfromGM0->data;
   //growfromWM->dim[4] = 2 ;  
   //nifti_update_dims_from_array(growfromGM);   // changing according sizes nt etc. 

    nifti_image * GMkoord0  = nifti_copy_nim_info(nim_input);
    nifti_image * GMkoord1  = nifti_copy_nim_info(nim_input);
    nifti_image * GMkoord2  = nifti_copy_nim_info(nim_input);
    nifti_image * GMkoord3  = nifti_copy_nim_info(nim_input);
    
    GMkoord0->datatype = NIFTI_TYPE_INT32; 
	GMkoord1->datatype = NIFTI_TYPE_INT32;
	GMkoord2->datatype = NIFTI_TYPE_INT32;
	GMkoord3->datatype = NIFTI_TYPE_INT32;

    GMkoord0->nbyper = sizeof(int);
	GMkoord1->nbyper = sizeof(int);
	GMkoord2->nbyper = sizeof(int);
	GMkoord3->nbyper = sizeof(int);

    GMkoord0->data = calloc(GMkoord0->nvox, GMkoord0->nbyper);
    GMkoord1->data = calloc(GMkoord1->nvox, GMkoord1->nbyper);
    GMkoord2->data = calloc(GMkoord2->nvox, GMkoord2->nbyper);
    GMkoord3->data = calloc(GMkoord3->nvox, GMkoord3->nbyper);

    int  *GMkoord0_data = (int *) GMkoord0->data;
    int  *GMkoord1_data = (int *) GMkoord1->data;
    int  *GMkoord2_data = (int *) GMkoord2->data;
    int  *GMkoord3_data = (int *) GMkoord3->data;


   
    nifti_image * equi_dist_layers  = nifti_copy_nim_info(nim_input);
	equi_dist_layers->datatype = NIFTI_TYPE_INT32;
	equi_dist_layers->nbyper = sizeof(int);
    equi_dist_layers->data = calloc(equi_dist_layers->nvox, equi_dist_layers->nbyper);
    int  *equi_dist_layers_data = (int *) equi_dist_layers->data;

   //nifti_image * equi_dist_layers = nifti_image_read(fin, 1);
   //short  *equi_dist_layers_data = (short *) equi_dist_layers->data;
   //equi_dist_layers->dim[4] = 1 ;  
   //nifti_update_dims_from_array(equi_dist_layers);   // changing according sizes nt etc.

//koordinaten


float x1g = 0.;
float y1g = 0.;
float x2g = 0.;
float y2g = 0.;
float x3g = 0.;
float y3g = 0.;

float dist (float x1, float y1, float x2, float y2) ; 
float angle (float a, float b, float c) ; 

cout << "bis hier 2 " << endl; 


// Reduce mask to contain only Areas close to the curface. 
cout << " select GM regions .... " << endl; 

int vinc = 700; // This is the distance from every voxel that the algorythm is applied on. Just to make it faster and not loop over all voxels.


float dist_i = 0.; 
float dist_min = 0.;
float dist_min1 = 0.;
float dist_min2 = 0.;
float dist_min3 = 0.;
float dist_max = 0.;
float dist_p1 = 0.;


int number_of_layers = 20 ; 


cout << " start growing  from WM .... " << endl; 



/// setting zero

    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
		 *(growfromWM0_data  + nxy*islice + nx*ix  + iy  ) = 0 ;
		 *(growfromWM1_data  + nxy*islice + nx*ix  + iy  ) = 0 ;
		 *(growfromGM0_data  + nxy*islice + nx*ix  + iy  ) = 0 ;
		 *(growfromGM1_data  + nxy*islice + nx*ix  + iy  ) = 0 ;
		 }
        }
      }
    



//////////////////////////////////
/////grow from  WM       /////////
//////////////////////////////////

int grow_vinc = 2 ; 



    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
	 	 if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 2 ) {
			*(growfromWM0_data  + nxy*islice + nx*ix  + iy  ) = 1.; 
			*(WMkoord0_data  + nxy*islice + nx*ix  + iy  ) = ix ; 
			*(WMkoord1_data  + nxy*islice + nx*ix  + iy  ) = iy ; 
		 }
        }
      }
    }

  for (int grow_i = 1 ; grow_i < vinc ; grow_i++ ){
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){

	dist_min2 = 10000.;
	  x1g = 0;
	  y1g = 0;
	   if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) ==  3  && *(growfromWM0_data  + nxy*islice + nx*ix  + iy  ) == 0 ){
	   	//cout << " true   " << *(growfromWM0_data  + nxy*islice + nx*ix  + iy  )<< endl; 
	    	for(int iy_i=max(0,iy-grow_vinc); iy_i<min(iy+grow_vinc,sizePhase); ++iy_i){
	     	 for(int ix_i=max(0,ix-grow_vinc); ix_i<min(ix+grow_vinc,sizeRead); ++ix_i){
			  if (*(growfromWM0_data  + nxy*islice + nx*ix_i  + iy_i  ) == (float)grow_i){
		 
			  dist_i = dist((float)ix,(float)iy,(float)ix_i,(float)iy_i); 
			  
			  if (dist_i < dist_min2 ){
			    dist_min2 = dist_i ; 
			    x1g = ix_i;
			    y1g = iy_i;
			    dist_p1 = dist_min2; 

			  }  
			}  
	  	 }
	  	}
		if ( dist_min2 < 1.4){
			//distDebug(0,islice,iy,ix) = dist_min2 ; 
 	    	*(growfromWM0_data  + nxy*islice + nx*ix  + iy  ) = (float)grow_i+1 ;
			*(WMkoord0_data  + nxy*islice + nx*ix  + iy  )  = *(WMkoord0_data  + nxy*islice + nx*(int)x1g  + (int)y1g  ); 
			*(WMkoord1_data  +  nxy*islice + nx*ix  + iy  ) = *(WMkoord1_data  + nxy*islice + nx*(int)x1g  + (int)y1g  );  
           
		}
		//cout << " ix   "  << ix << " iy   "  << iy  << "    " << *(WMkoord0_data  + nxy*islice + nx*(int)x1g  + (int)y1g  )<< endl; 
	   }
	   

        }
      }
    }
 }


//////////////////////////////////
/////grow from  CSF       /////////
//////////////////////////////////

cout << " start growing from CSF .... " << endl; 

    
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
	 	if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 1 ) {
			*(growfromGM0_data  + nxy*islice + nx*ix  + iy  ) = 1.; 
			*(GMkoord0_data  + nxy*islice + nx*ix  + iy  ) = ix ; 
			*(GMkoord1_data  + nxy*islice + nx*ix  + iy  ) = iy ; 
		}
        }
      }
    }    
    
    

  for (int grow_i = 1 ; grow_i < vinc ; grow_i++ ){
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){

	dist_min2 = 10000.;
	  x1g = 0.;
	  y1g = 0;
	   if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 3  && *(growfromGM0_data  + nxy*islice + nx*ix  + iy  )  == 0 ){
	    for(int iy_i=max(0,iy-grow_vinc); iy_i<min(iy+grow_vinc,sizePhase); ++iy_i){
	      for(int ix_i=max(0,ix-grow_vinc); ix_i<min(ix+grow_vinc,sizeRead); ++ix_i){
		if (*(growfromGM0_data  + nxy*islice + nx*ix_i  + iy_i  )  == (float)grow_i){
		 
		  dist_i = dist((float)ix,(float)iy,(float)ix_i,(float)iy_i); 
		  if (dist_i < dist_min2 ){
		    dist_min2 = dist_i ; 
		    x1g = ix_i;
		    y1g = iy_i;
		    dist_p1 = dist_min2; 
		  }  
		}  
	      }
	    }
		if ( dist_min2 < 1.4){
 	    	*(growfromGM0_data  + nxy*islice + nx*ix  + iy  )  = (float)grow_i+1 ;

		    *(GMkoord0_data  + nxy*islice + nx*ix  + iy  )         = *(GMkoord0_data  +          nxy*islice + nx*(int)x1g  + (int)y1g  ); 
			*(GMkoord1_data  + nxy*islice + nx*ix  + iy  )         = *(GMkoord1_data  +          nxy*islice + nx*(int)x1g  + (int)y1g  );  
		}
	   }

        }
      }
    }
 }


/////////////////////////////////////////////////////////////////////////////////////////////////////
///// wabble accross neigbouring voexles of closest WM to account for Pytagoras errors      /////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

cout << " correct for pytagoras error .... " << endl; 


    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
			*(growfromWM1_data +   nxy*islice + nx*ix  + iy  ) =   *(growfromWM0_data  +  nxy*islice + nx*ix  + iy  ) ; 
			*(WMkoord2_data  +  nxy*islice + nx*ix  + iy  ) = *(WMkoord0_data  +  nxy*islice + nx*ix  + iy  )  ; 
			*(WMkoord3_data  +  nxy*islice + nx*ix  + iy  ) = *(WMkoord1_data  +  nxy*islice + nx*ix  + iy  )  ; 

        }
      }
    }

  for (int grow_i = 1 ; grow_i < vinc ; grow_i++ ){
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
	   if (*(WMkoord1_data   +  nxy*islice + nx*ix  + iy  )  != 0 ){
		dist_min2 = 10000.;
	  	x1g = 0;
	  	y1g = 0;

	    	for(int iy_i=max(0,(*(WMkoord3_data   +  nxy*islice + nx*ix  + iy  ))-grow_vinc);   iy_i<min((*(WMkoord3_data  +  nxy*islice + nx*ix  + iy  ))+grow_vinc,sizePhase); ++iy_i){
	    	  for(int ix_i=max(0,(*(WMkoord2_data   +  nxy*islice + nx*ix  + iy  ))-grow_vinc); ix_i<min((*(WMkoord2_data  +  nxy*islice + nx*ix  + iy  ))+grow_vinc,sizeRead); ++ix_i){
			    if (*(nim_input_data  + nxy*islice + nx*ix_i  + iy_i  )  == 2){
			 
			  dist_i =  dist((float)ix,(float)iy,(float)ix_i,(float)iy_i);
			  if (dist_i < dist_min2 ){
			    dist_min2 = dist_i ; 
			    x1g = ix_i;
			    y1g = iy_i;
			    dist_p1 = dist_min2; 
			  }  
			}  
		   }
		}
 	    *(growfromWM1_data  + nxy*islice + nx*ix  + iy  ) =  dist((float)ix,(float)iy,(float)x1g,(float)y1g); 
		*(WMkoord2_data   +  nxy*islice + nx*ix  + iy  ) = *(WMkoord2_data  +   nxy*islice + nx*(int)x1g  + (int)y1g  ); 
		*(WMkoord3_data   +  nxy*islice + nx*ix  + iy  ) = *(WMkoord3_data  +   nxy*islice + nx*(int)x1g  + (int)y1g  );  
	   }
        }
      }
    }
 }

cout << " runing until here .... " << endl; 


/////////////////////////////////////////////////////////////////////////////////////////////////////
///// wabble accross neigbouring voexles of closest GM to account for Pytagoras errors      /////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

    
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
			*(growfromGM1_data  +  nxy*islice + nx*ix  + iy  )  =   *(growfromGM0_data        +  nxy*islice + nx*ix  + iy  ) ; 
			*(GMkoord2_data      +  nxy*islice + nx*ix  + iy  ) =   *(GMkoord0_data  +  nxy*islice + nx*ix  + iy  )  ; 
			*(GMkoord3_data    +  nxy*islice + nx*ix  + iy  )   =   *(GMkoord1_data    +  nxy*islice + nx*ix  + iy  )  ; 

        }
      }
    }

cout << " runing also until here .... " << endl; 




  for (int grow_i = 1 ; grow_i < vinc ; grow_i++ ){
    for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
	   if (*(GMkoord1_data  +  nxy*islice + nx*ix  + iy  )  != 0 ){
		dist_min2 = 10000.;
	  	x1g = 0;
	  	y1g = 0;

	    	for(int iy_i=max(0,(*(GMkoord3_data    +  nxy*islice + nx*ix  + iy  ))-grow_vinc);   iy_i<min((*(GMkoord3_data  +  nxy*islice + nx*ix  + iy  ))+grow_vinc,sizePhase); ++iy_i){
	    	  for(int ix_i=max(0,(*(GMkoord2_data  +  nxy*islice + nx*ix  + iy  ))-grow_vinc);   ix_i<min((*(GMkoord2_data  +  nxy*islice + nx*ix  + iy  ))+grow_vinc,sizeRead);  ++ix_i){
			   if (*(nim_input_data  + nxy*islice + nx*ix_i  + iy_i  )  == 1){
			 
			  dist_i =  dist((float)ix,(float)iy,(float)ix_i,(float)iy_i);
			  if (dist_i < dist_min2 ){
			    dist_min2 = dist_i ; 
			    x1g = ix_i;
			    y1g = iy_i;
			    dist_p1 = dist_min2; 
			  }  
			}  
		   }
		}
 	    *(growfromGM1_data  + nxy*islice + nx*ix  + iy  ) =  dist((float)ix,(float)iy,(float)x1g,(float)y1g); 
		*(GMkoord2_data   +  nxy*islice + nx*ix  + iy  ) = *(GMkoord2_data  +   nxy*islice + nx*(int)x1g  + (int)y1g  ); 
		*(GMkoord3_data   +  nxy*islice + nx*ix  + iy  ) = *(GMkoord3_data   +   nxy*islice + nx*(int)x1g  + (int)y1g  );  
	   }
        }
      }
    }
 }

cout << " runing also until here 3 .... " << endl; 

int   GMK2_i, GMK3_i, WMK2_i, WMK3_i ; 
float GMK2_f, GMK3_f, WMK2_f, WMK3_f, ix_f, iy_f ; 

for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
		if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 3){
	  		//equi_dist_layers(0,islice,iy,ix)                 =  19 * (1- dist((float)ix,(float)iy,(float)GMkoord(2,islice,iy,ix)                     ,(float)GMkoord(3,islice,iy,ix)                            )/ (dist((float)ix,(float)iy,(float)GMkoord(2,islice,iy,ix)                          ,(float)GMkoord(3,islice,iy,ix)                            ) + dist((float)ix,(float)iy,(float)WMkoord(2,islice,iy,ix)                              ,(float)WMkoord(3,islice,iy,ix)                        ) )) + 2  ;
            
            GMK2_i  = *(GMkoord2_data + nxy*islice + nx*ix + iy) ; 
            GMK3_i  = *(GMkoord3_data + nxy*islice + nx*ix + iy) ;       
            WMK2_i  = *(WMkoord2_data + nxy*islice + nx*ix + iy) ; 
            WMK3_i  = *(WMkoord3_data + nxy*islice + nx*ix + iy) ;    
            GMK2_f  = (float)GMK2_i; 
            GMK3_f  = (float)GMK3_i;      
            WMK2_f  = (float)WMK2_i; 
            WMK3_f  = (float)WMK3_i;   
            ix_f    = (float)ix;  
            iy_f    = (float)iy;  

            // cout << " rix_f,iy_f,GMK2_f,GMK3_f " <<  "   "  <<  ix_f <<  "   "  <<iy_f <<  "   "  <<GMK2_f <<  "   "  <<*(GMkoord2_data + nxy*islice + nx*ix + iy)<< endl; 
            *(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) =  19 * (1- dist((float)ix,(float)iy,GMK2_f,GMK3_f )/ (dist((float)ix,(float)iy,GMK2_f,GMK3_f) + dist((float)ix,(float)iy,WMK2_f,WMK3_f) )) + 2  ;
            //*(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) =  100 * dist(ix_f,iy_f,GMK2_f,GMK3_f ) ;

            	}
        }
      }
    }

 cout << " runing also until here  4.... " << endl; 

// Cleaning negative layers and layers ov more than 20

for(int islice=0; islice<sizeSlice; ++islice){  
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead-0; ++ix){
		if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 1 && *(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) == 0  ){
	  		*(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) = 21 ;
            	}
		if (*(nim_input_data  + nxy*islice + nx*ix  + iy  ) == 2 && *(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) == 0  ){
	  		*(equi_dist_layers_data  +  nxy*islice + nx*ix  + iy  ) = 1 ;
            	}
        }
      }
    }
 cout << " runing also until here  4.5... " << endl; 

//equi_dist_layers.autowrite("equi_dist_layers.nii", wopts, &prot);

 cout << " runing also until here  5.... " << endl; 


     // output file name       
  const char  *fout_4="equi_dist_layers.nii" ;
  if( nifti_set_filenames(equi_dist_layers, fout_4 , 1, 1) ) return 1;
  nifti_image_write( equi_dist_layers );

 // const char  *fout_5="debug_ing.nii" ;
 // if( nifti_set_filenames(growfromWM0, fout_5 , 1, 1) ) return 1;
 // nifti_image_write( growfromWM0 );
  
 // const char  *fout_6="kootrGM.nii" ;
 // if( nifti_set_filenames(GMkoord2, fout_6 , 1, 1) ) return 1;
 // nifti_image_write( GMkoord2 );

 // koord.autowrite("koordinaten.nii", wopts, &prot);
  return 0;
}



  float dist (float x1, float y1, float x2, float y2) {
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  }

  float angle (float a, float b, float c) {
	if (a*a+b*b-c*c <= 0 ) return 3.141592 ;
    	else return acos((a*a+b*b-c*c)/(2.*a*b));
  }

