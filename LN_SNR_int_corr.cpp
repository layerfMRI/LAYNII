#include <stdio.h>
#include "nifti2_io.h"
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
      "       -Nulled_file  INFILE      : specify input dataset\n"
      "       -Not_Nulled_file  INFILE      : specify input dataset\n"
      "       -output OUTFILE     : specify output dataset\n"
      "       -verb LEVEL         : set the verbose level to LEVEL\n"
      "       -cutoff    value    : set a cutof\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{
   nifti_image * nim_n=NULL, * nim_nn=NULL;
   char        * fin_n=NULL, * fin_nn=NULL; 
   int         cutoff,  ac, disp_float_eg=0;

   if( argc < 2 ) return show_help();   /* typing '-help' is sooo much work */

   /* process user options: 4 are valid presently */
   for( ac = 1; ac < argc; ac++ ) {
      if( ! strncmp(argv[ac], "-h", 2) ) {
         return show_help();
      }
      else if( ! strcmp(argv[ac], "-disp_float_example") ) {
         disp_float_eg = 1;
      }
      else if( ! strcmp(argv[ac], "-Nulled_file") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -Nulled_file\n");
            return 1;
         }
         fin_n = argv[ac];  /* no string copy, just pointer assignment */
      }
      else if( ! strcmp(argv[ac], "-Not_Nulled_file") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -Nulled_file\n");
            return 1;
         }
         fin_nn = argv[ac];  /* no string copy, just pointer assignment */
      }
      else if( ! strcmp(argv[ac], "-cutoff") ) {
         if( ++ac >= argc ) {
            fprintf(stderr, "** missing argument for -input\n");
            return 1;
         }
         cutoff = atoi(argv[ac]);  /* no string copy, just pointer assignment */
         cout << " cutoff is " << cutoff << endl;
      }
      else {
         fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
         return 1;
      }
   }

   if( !fin_n  ) { fprintf(stderr, "** missing option '-Nulled'\n");  return 1; }
   /* read input dataset, including data */
   nim_n = nifti_image_read(fin_n, 1);
   if( !nim_n ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin_n);
      return 2;
   }
   
   if( !fin_nn  ) { fprintf(stderr, "** missing option '-Not_Nulled'\n");  return 1; }
   /* read input dataset, including data */
   nim_nn = nifti_image_read(fin_nn, 1);
   if( !nim_nn ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin_nn);
      return 2;
   }
   

   // get dimsions of input 
   int sizeSlice = nim_n->nz ; 
   int sizePhase = nim_n->nx ; 
   int sizeRead = nim_n->ny ; 
   int nrep =  nim_n->nt; 
   
   cout << sizeSlice << " slices    " <<  sizePhase << " PhaseSteps     " <<  sizeRead << " Read steps    " <<  nrep << " timesteps "  << endl; 

   int nx =  nim_n->nx;
   int nxy = nim_n->nx * nim_n->ny;
   int nxyz = nim_n->nx * nim_n->ny * nim_n->nz;



  // allocating an additional nii 
   float  *nim_nn_data = (float *) nim_nn->data;
   float  *nim_n_data = (float *) nim_n->data;
   cout << " bis hier2   " << *(nim_n_data + nxyz*2 + nxy*8 + nx*64  + 64  )<<  endl; 

   nifti_image * tSNR_abs_nulled = nim_n; 
   float  *tSNR_abs_nulled_data = (float *) tSNR_abs_nulled->data;
   tSNR_abs_nulled->dim[4] = 1 ;  
   nifti_update_dims_from_array(tSNR_abs_nulled);   // changing according sizes nt etc. 
   
   nifti_image * tSNR_nulled = nim_n; 
   float  *tSNR_nulled_data = (float *)tSNR_nulled->data;
   tSNR_nulled->dim[4] = 1 ;  
   nifti_update_dims_from_array(tSNR_nulled);   // changing according sizes nt etc. 
   
   nifti_image * mean = nim_n; 
   float  *mean_data = (float *)mean->data;
   mean->dim[4] = 1 ;  
   nifti_update_dims_from_array(mean);   // changing according sizes nt etc. 

int N = nrep/2 ; 
double vec_n[N]  ;
double vec_nn[N]  ;


cout << " nrep " <<  nrep  << endl; 




    for(int islice=0; islice<sizeSlice; ++islice){
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
	  // if ( *(nim_n_data + nxyz*1 + nxy*islice + nx*ix  + iy  )  > cutoff &&  *(nim_n_data + nxyz*1 + nxy*islice + nx*ix  + iy  )  < 100000000){
		for(int timestep=0; timestep<N-40 ; timestep = timestep + 2 ) {
		

		
			vec_n[(int)timestep]  = *(nim_n_data + nxyz*(timestep+1) + nxy*islice + nx*ix  + iy  ); 
			vec_n[(int)timestep+1]  = *(nim_nn_data + nxyz*(timestep+2) + nxy*islice + nx*ix  + iy  );
			vec_nn[(int)timestep/2] = *(nim_n_data + nxyz*(timestep+1) + nxy*islice + nx*ix  + iy  ); 
			//cout << " islice " <<  islice << " sizePhase " <<  iy << " sizeRead " <<  ix << " timestep " <<  timestep << "   " << *(nim_n_data + nxyz*(timestep+1) + nxy*islice + nx*ix  + iy  )<<  endl; 

	   	 }
        *(tSNR_nulled_data + nxyz*0 + nxy*islice + nx*ix  + iy  ) = (int) gsl_stats_mean (vec_n, 1, N-2) /gsl_stats_sd_m(vec_n, 1, N-2,  gsl_stats_mean (vec_n, 1, N-2));
      //  *(tSNR_abs_nulled  + nxyz*0 + nxy*islice + nx*ix  + iy  ) = (float) gsl_stats_mean (vec_nn, 1, (N-2)/2) /gsl_stats_sd_m(vec_nn, 1, (N-2)/2,  gsl_stats_mean (vec_nn, 1, (N-2)/2));
	   //tSNR_nulled(0,islice,iy,ix) = gsl_stats_sd_m(vec_n, 1, N-2,  gsl_stats_mean (vec_n, 1, N-2));
		//cout << " tSNR_nulled(0,islice,iy,ix) " << tSNR_nulled(0,islice,iy,ix) << endl; 
		//gsl_vector_set_zero (my_gsl_vec_n);
		
				//cout << " islice " <<  islice << " sizePhase " <<  iy << " sizeRead " <<  ix << " timestep " <<  timestep << endl; 

	    *(mean_data+ nxyz*0 + nxy*islice + nx*ix  + iy  ) = (int) gsl_stats_mean (vec_n, 1, N-2) ;
     //     }
	// else 	{
	//	*(tSNR_nulled_data + nxyz*0 + nxy*islice + nx*ix  + iy  ) = 0.; 
	//	}
        }
      } 
     }


     
     
     
/*
// Intensity correction: 


  Data<float,4> tSNR_intens;
  tSNR_intens.resize(1,sizeSlice,sizePhase,sizeRead);
  tSNR_intens=0.0;

  
  Data<float,4> mean_intens;
  mean_intens.resize(1,sizeSlice,sizePhase,sizeRead);
  mean_intens=0.0;


double N_voxels_sl_tsnr = 0.; 
double intens_sl_tsnr = 0.; 
double intens_sl_mean = 0.; 

    for(int islice=0; islice<sizeSlice; ++islice){
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
	    
	  if (mean(0,islice,iy,ix) > cutoff){
	    
	    intens_sl_tsnr = intens_sl_tsnr + tSNR_nulled(0,islice,iy,ix) ; 
	    intens_sl_mean = intens_sl_mean + mean(0,islice,iy,ix) ; 
	    N_voxels_sl_tsnr  = N_voxels_sl_tsnr + 1; 
	  
	  }
	}
      }
      
      cout << "slice " << islice << " has "<<   N_voxels_sl_tsnr<< " voxels  with mean  intensity " << intens_sl_tsnr/N_voxels_sl_tsnr << endl; 
      
      for(int iy=0; iy<sizePhase; ++iy){
        for(int ix=0; ix<sizeRead; ++ix){
	    
	  if (mean(0,islice,iy,ix) > cutoff){
	    
	   tSNR_intens(0,islice,iy,ix)  =  tSNR_nulled(0,islice,iy,ix) / (intens_sl_tsnr / N_voxels_sl_tsnr); 
	   mean_intens(0,islice,iy,ix)  =  mean(0,islice,iy,ix)        / (intens_sl_mean / N_voxels_sl_tsnr); 

	  
	  }
	}
      }
      
         N_voxels_sl_tsnr = 0.; 
	 intens_sl_tsnr = 0.; 
	 intens_sl_mean = 0.;  
    }    
*/

//cout << " bis hier3 " << endl; 

  const char  *fout_3="tSNR_abs.nii" ;
  if( nifti_set_filenames(tSNR_abs_nulled, fout_3 , 1, 1) ) return 1;
  nifti_image_write( tSNR_abs_nulled );
  
  const char  *fout_4="tSNR_nulled.nii" ;
  if( nifti_set_filenames(tSNR_nulled, fout_4 , 1, 1) ) return 1;
  nifti_image_write( tSNR_nulled );

  const char  *fout_5="mean.nii" ;
  if( nifti_set_filenames(mean, fout_5 , 1, 1) ) return 1;
  nifti_image_write( mean );


  const char  *fout_6="nulled_as_in.nii" ;
  if( nifti_set_filenames(nim_n, fout_6 , 1, 1) ) return 1;
  nifti_image_write( nim_n );
  // tSNR_nulled.autowrite("T1_instability_"+filename1, wopts, &prot);
  // tSNR_abs_nulled.autowrite("tSNR_"+filename1, wopts, &prot);
  // mean.autowrite("MEAN_"+filename1, wopts, &prot);
//   tSNR_intens.autowrite("T1_instability_intens"+filename1, wopts, &prot);
//   mean_intens.autowrite("MEAN_intens"+filename1, wopts, &prot);

//cout << " bis hier4 " << endl; 

  return 0;

}
