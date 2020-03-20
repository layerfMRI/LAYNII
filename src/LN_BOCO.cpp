
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_BOCO: This program does BOLD correction in SS-SI VASO. It does \n"
    "         the division of nulled and not nulled imaged. \n"
    "\n"
    "Usage:\n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii \n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -shift \n"
    "    LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 24 \n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -Nulled    : Nulled (VASO) time series that needs to be BOLD \n"
    "               : corrected.\n"
    "    -BOLD      : Reference BOLD time series without a VASO contrast.\n"
    "    -shift     : (Optional) Estimate the correlation of BOLD and VASO \n"
    "                 for temporal shifts.\n"
    "    -trialBOCO : First average trials and then do the BOLD correction. \n"
    "                 The parameter is the trial duration in TRs.\n"
    "\n"
    "Notes:\n"
    "    - Here it is assumed that BOLD and VASO refer to the double TR: \n"
    "        3dUpsample -overwrite -datum short -prefix Nulled_intemp.nii -n 2 -input Nulled.nii \n"
    "        3dUpsample -overwrite -datum short -prefix BOLD_intemp.nii -n 2 -input BOLD.nii \n"
    "    - Here I assume that they have the same spatiotemporal dimensions. \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char* fin_1 = NULL, * fin_2 = NULL;
    int ac, shift = 0;
    int trialdur = 0;
    if (argc < 2) return show_help();

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-Nulled")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Nulled\n");
                return 1;
            }
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-BOLD")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -BOLD\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-trialBOCO")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -trialBOCO\n");
                return 1;
            }
            trialdur = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-shift")) {
            shift = 1;
            cout << "Do a correlation analysis with temporal shifts."  << endl;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-Nulled'.\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-BOLD'.\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'.\n", fin_2);
        return 2;
    }

    log_welcome("LN_BOCO");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int size_time = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;
    const int nr_voxels = size_time * size_z * size_y * size_x;

    // ========================================================================
    // Fix datatype issues
    nifti_image *nii_nulled = copy_nifti_as_float32(nii1);
    float *nii_nulled_data = static_cast<float*>(nii_nulled->data);
    nifti_image *nii_bold = copy_nifti_as_float32(nii2);
    float *nii_bold_data = static_cast<float*>(nii_bold->data);

    // Allocate new nifti
    nifti_image *nii_boco_vaso = copy_nifti_as_float32(nii1);
    float *nii_boco_vaso_data = static_cast<float*>(nii_boco_vaso->data);

    // ========================================================================
    // AVERAGE across Trials
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_boco_vaso_data + i) = *(nii_nulled_data + i)
                                    / (*(nii_bold_data + i));
    }

    // Clean VASO values that are unrealistic
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_boco_vaso_data + i) <= 0) {
            *(nii_boco_vaso_data + i) = 0;
        }
        if (*(nii_boco_vaso_data + i) >= 5) {
            *(nii_boco_vaso_data + i) = 5;
        }
    }

    // ========================================================================
    // Shift
    if (shift == 1) {
        nifti_image* correl_file  = nifti_copy_nim_info(nii_nulled);
        correl_file->nt = 7;
        correl_file->nvox = nii_nulled->nvox / size_time *7;
        correl_file->datatype = NIFTI_TYPE_FLOAT32;
        correl_file->nbyper = sizeof(float);
        correl_file->data = calloc(correl_file->nvox, correl_file->nbyper);
        float* correl_file_data = static_cast<float*>(correl_file->data);

        double vec_file1[size_time];
        double vec_file2[size_time];

        for (int shift = -3; shift <= 3; ++shift) {
            cout << "  Calculating shift = " << shift << endl;
            for (int j = 0; j != size_z * size_y * size_x; ++j) {
                for (int t = 3; t < size_time-3; ++t) {
                    *(nii_boco_vaso_data + nxyz * t + j) =
                        *(nii_nulled_data + nxyz * t + j)
                        / *(nii_bold_data + nxyz * (t + shift) + j);
                }
                for (int t = 0; t < size_time; ++t) {
                    vec_file1[t] = *(nii_boco_vaso_data + nxyz * t + j);
                    vec_file2[t] = *(nii_bold_data + nxyz * t + j);
                }
                *(correl_file_data + nxyz * (shift + 3) + j) =
                    ren_correl(vec_file1, vec_file2, size_time);
            }
        }

        // Get back to default
        for (int i = 0; i != nr_voxels; ++i) {
            *(nii_boco_vaso_data + i) = *(nii_nulled_data + i)
                                        / *(nii_bold_data + i);
        }

        // Clean VASO values that are unrealistic
        for (int i = 0; i != nr_voxels; ++i) {
            if (*(nii_boco_vaso_data + i) <= 0) {
                *(nii_boco_vaso_data + i) = 0;
            }
            if (*(nii_boco_vaso_data + i) >= 2) {
                *(nii_boco_vaso_data + i) = 2;
            }
        }
        save_output_nifti(fin_1, "correlated", correl_file, false);
    }

    // ========================================================================
    // Trial average
if (trialdur!=0) {

cout << " I will also do the BOLD correction after the trial average " << endl; 

   
   
   int nrep = size_time ; 
   int sizeRead= size_y ;
   int sizePhase = size_x ;
   int sizeSlice = size_z ;
   
   cout << " Trial duration is " <<trialdur << " this means there are " << (float)nrep/(float)trialdur <<  " trials recorted here " << endl; 


   int numberofTrials = nrep/trialdur ; 
   // Trial averave file
    nifti_image * triav_file    = nifti_copy_nim_info(nii1);
    triav_file->nt 				= trialdur	; 
    triav_file->nvox 			= nii1->nvox / nrep * trialdur; 
    triav_file->datatype 		= NIFTI_TYPE_FLOAT32; 
    triav_file->nbyper 			= sizeof(float);
    triav_file->data 			= calloc(triav_file->nvox, triav_file->nbyper);
    float  *nii_triavg_data 	= (float *) triav_file->data;
    
    
    nifti_image * triav_B_file  = nifti_copy_nim_info(nii1);
    triav_B_file->nt 			= trialdur	; 
    triav_B_file->nvox 			= nii1->nvox / nrep * trialdur; 
    triav_B_file->datatype 		= NIFTI_TYPE_FLOAT32; 
    triav_B_file->nbyper 		= sizeof(float);
    triav_B_file->data 			= calloc(triav_B_file->nvox, triav_B_file->nbyper);
    float  *nii_triavg_B_data 	= (float *) triav_B_file->data;
    
    float AV_Nulled[trialdur] ;
    float AV_BOLD[trialdur]   ; 


    
    
	  for(int islice=0; islice<sizeSlice; ++islice){  
	      for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
	              for(int it=0; it<trialdur; ++it){  
					AV_Nulled[it] = 0 ;
					AV_BOLD  [it] = 0 ;
				  }  
	           	  for(int it=0; it<trialdur*numberofTrials; ++it){  
        		  		AV_Nulled[it%trialdur] = AV_Nulled[it%trialdur] + (*(nii_nulled_data  + nxyz *(it) +  nxy*islice + nx*ix  + iy  ))/numberofTrials ;	
        		  		AV_BOLD[it%trialdur]   = AV_BOLD[it%trialdur]   + (*(nii_bold_data  + nxyz *(it) +  nxy*islice + nx*ix  + iy  ))/numberofTrials ;	
			      }         
			      
			      for(int it=0; it<trialdur; ++it){  
			        *(nii_triavg_data    + nxyz *it +  nxy*islice + nx*ix  + iy  ) = AV_Nulled[it]/AV_BOLD[it] ;
			        *(nii_triavg_B_data  + nxyz *it +  nxy*islice + nx*ix  + iy  ) = AV_BOLD[it] ;

				  } 
           } 
	    }
	  }
	
     // clean VASO values that are unrealistic
    for(int islice=0; islice<sizeSlice; ++islice){  
	   for(int iy=0; iy<sizePhase; ++iy){
	        for(int ix=0; ix<sizeRead; ++ix){
	        	for(int it=0;  it<trialdur; ++it){  
	
	            	if (*(nii_triavg_data + nxyz*it + nxy*islice + nx*ix  + iy  ) <= 0) {
	            		*(nii_triavg_data + nxyz*it + nxy*islice + nx*ix  + iy  ) = 0 ;
	            		}
	            	if (*(nii_triavg_data + nxyz*it + nxy*islice + nx*ix  + iy  ) >= 2) {
	            		*(nii_triavg_data + nxyz*it + nxy*islice + nx*ix  + iy  ) = 2 ;
	            		}

	            }
        	}	
	    }
	  }


  const char  *fout_trial="VASO_trialAV_LN.nii" ;
  if( nifti_set_filenames(triav_file, fout_trial , 1, 1) ) return 1;
  nifti_image_write( triav_file );

  const char  *fout_trial_BOLD="BOLD_trialAV_LN.nii" ;
  if( nifti_set_filenames(triav_B_file, fout_trial_BOLD , 1, 1) ) return 1;
  nifti_image_write( triav_B_file );


}// Trial Average loop closed
    
    const char  *fout_5="VASO_LN.nii" ;
    if( nifti_set_filenames(nii_boco_vaso, fout_5 , 1, 1) ) return 1;
    nifti_image_write( nii_boco_vaso );
    
//  RENZO needs to include an output option here. For the time beeing, and for the 
//    save_output_nifti(fin_2, "VASO", nii_boco_vaso, true);

    cout << "  Finished." << endl;
    return 0;
}
