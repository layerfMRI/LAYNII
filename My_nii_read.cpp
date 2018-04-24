//compile with with "make My_nii_read"
//execute with ./My_nii_read -input input_example.nii -output output.nii -cutoff 3

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

int show_help( void )
{
   printf(
      "My_nii_read: short exmample of reading/writing NIfTI2\n"
      "\n"
      "    This program is to demonstrate how to read a NIfTI-2 dataset,\n"
      "    set output filenames and write a NIfTI-2 dataset, all via the\n"
      "    standard NIfTI C library.\n"
      "\n"
      "    basic usage: ./My_nii_read -input input_example.nii -output output.nii -cutoff 3\n"
      "\n"
      "    options:\n"
      "\n"
      "       -help               : show this help\n"
      "       -disp_float_example : show some voxel's data\n"
      "       -input  INFILE      : specify input dataset\n"
      "       -output OUTFILE     : specify output dataset\n"
      "       -verb LEVEL         : set the verbose level to LEVEL\n"
      "       -cutoff    value    : set a cutof\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{

   nifti_image * nim_input=NULL;
   char        * fin=NULL, * fout=NULL;
   int         cutoff,  ac, disp_float_eg=0;

   if( argc < 2 ) return show_help();   // typing '-help' is sooo much work 

   // process user options: 4 are valid presently 
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-disp_float_example") ) {
         disp_float_eg = 1;
      }
      else if( ! strcmp(argv[ac], "-input") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         fin = argv[ac];  // no string copy, just pointer assignment 
      }
      else if( ! strcmp(argv[ac], "-output") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -output\n");
            return 2;
         }
         fout = argv[ac];
      }
      else if( ! strcmp(argv[ac], "-verb") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -verb\n");
            return 2;
         }
         nifti_set_debug_level(atoi(argv[ac]));
      }
      else if( ! strcmp(argv[ac], "-cutoff") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         cutoff = atoi(argv[ac]);  // no string copy, just pointer assignment 
         cout << " cutoff is " << cutoff << endl;
      }
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !fin  ) { fprintf(stderr, "** missing option '-input'\n");  return 1; }
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


   if( !fout ) { fprintf(stderr, "-- no output requested \n"); return 0; }
   // assign nifti_image fname/iname pair, based on output filename
   //   (request to 'check' image and 'set_byte_order' here) 
   if( nifti_set_filenames(nim_input, fout, 1, 1) ) return 1;

	// get access to data of nim_input 
    float  *nim_input_data = (float *) nim_input->data;
    
    
   // allocating an additional nii 
   nifti_image * nim_output = nifti_image_read(fin, 1);
   //nifti_image * nim_output =nim_input; 
   float  *nim_output_data = (float *) nim_output->data;

   // Changing the number of TRs. So far it only works to reduce the number of TRs. when its larger, it doesn't work 
   nim_output->dim[4] = 1 ;  
   nifti_update_dims_from_array(nim_output);   // changing according sizes nt etc. 
   //cout << " has valid dims =  " <<  nifti_nim_has_valid_dims    (nim_1, 0) << endl;
   //cout << " valid_nifti_extensions(nim_1) " <<  valid_nifti_extensions(nim_1) << endl;	 

// get input as 4D array
    cout <<  " loading to array works0 "  << endl; 

    cout <<  " loading to array works1 "  << endl; 

    cout <<  " loading to array works 2"  << endl; 


 //My individual analysis 
 cout << " ok"  << endl; 

  cout << " not ok"  << endl; 


 
  for(int timestep=0; timestep<nrep; ++timestep){
    for(int islice=0; islice<sizeSlice; ++islice){
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
        
        	if ( *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix  + iy  )  > cutoff )  {
        		*(nim_output_data + nxyz*timestep + nxy*islice + nx*ix  + iy  ) =   *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix  + iy  )  ;  
        	}  
        	else {
        		*(nim_output_data + nxyz*timestep + nxy*islice + nx*ix  + iy  ) = -0.9 *   *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix  + iy  )  ;  
        	}
        }  
      } 
    }
  }

    

   // output file name       
  const char  *fout_3="output_manipulated.nii" ;
  if( nifti_set_filenames(nim_output, fout_3 , 1, 1) ) return 1;
  nifti_image_write( nim_output );
  
   // writing out input file 
   // if we get here, write the output dataset 
   
  // if( nifti_set_filenames(nim_input, fout , 1, 1) ) return 1;
   nifti_image_write( nim_input );
  //  and clean up memory 
   nifti_image_free( nim_output );
    //nifti_image_free( nim_input );



return 0;

}



