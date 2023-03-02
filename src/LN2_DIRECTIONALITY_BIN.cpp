#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN2_DIRECTIONALITY_BIN: Estimates measure of columnarity \n"
    "                        and laminarity for binarizes maps \n"
    "                        with binarizes layers and columns as input.\n"
    "\n"
    "Usage:\n"
    "    LN2_DIRECTIONALITY_BIN -input binarized_map.nii -layers layers.nii.gz -columns columns.nii.gz \n"
    "\n"
    "test usage in the test_data folder: \n"
    "    ../LN2_DIRECTIONALITY_BIN -input binarized_map.nii -layers layers.nii.gz -columns columns.nii.gz \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -input       : Nifti (.nii) binarized map\n"
    "    -kernel_size : (Optional) Use an odd positive integer (default 11).\n"
    "    -output      : (Optional) Output filename, including .nii or\n"
    "                   .nii.gz, and path if needed. Overwrites existing files.\n"
    "                   If not given, the prefix 'fPSF' is added.\n"
    "\n"
    "Notes:\n"
    "    This is written foir Richard as a side project.  \n"
    "\n");
    return 0;
}




int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ;
    char *fin = NULL;
    char *fin_layers = NULL, *fin_columns = NULL;
    int ac;
    int kernel_size = 11; // This is the maximal number of layers. I don't know how to allocate it dynamically. this should be an odd number. That is smaller than half of the shortest matrix size to make sense
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-layers")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -layers\n");
                return 1;
            }
            fin_layers = argv[ac];
        } else if (!strcmp(argv[ac], "-columns")) {
            if (++ac >= argc) {
                fprintf(stderr, " ** missing argument for -columns\n");
                return 1;
            }
            fin_columns = argv[ac];
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        }
    }
    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_columns) {
        fprintf(stderr, " ** missing option '-columns'\n");
        return 1;
    }
    if (!fin_layers) {
        fprintf(stderr, " ** missing option '-layers'\n");
        return 1;
    }



    // Read input dataset
    nifti_image* nim_column_r = nifti_image_read(fin_columns, 1);
    if (!nim_column_r) {
        fprintf(stderr, " ** failed to read NIfTI from '%s'\n", fin_columns);
        return 2;
    }
    nifti_image* nim_layers_r = nifti_image_read(fin_layers, 1);
    if (!nim_layers_r) {
        fprintf(stderr, " ** failed to read NIfTI from '%s'\n", fin_layers);
        return 2;
    }
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN2_DIRECTIONALITY_BIN");
    log_nifti_descriptives(nii_input);
    log_nifti_descriptives(nim_layers_r);
    log_nifti_descriptives(nim_column_r);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    nifti_image* nim_layers = copy_nifti_as_int32(nim_layers_r);
    int32_t* nim_layers_data = static_cast<int32_t*>(nim_layers->data);

    nifti_image* nim_columns = copy_nifti_as_int32(nim_column_r);
    int32_t* nim_columns_data = static_cast<int32_t*>(nim_columns->data);

    // Allocate new nifti images
    nifti_image * nii_laminarity = nifti_copy_nim_info(nii);
    nii_laminarity->nt = 1;
    nii_laminarity->nvox = nii->nvox / size_time;
    nii_laminarity->datatype = NIFTI_TYPE_FLOAT32;
    nii_laminarity->nbyper = sizeof(float);
    nii_laminarity->data = calloc(nii_laminarity->nvox, nii_laminarity->nbyper);
    float* nii_laminarity_data = static_cast<float*>(nii_laminarity->data);
    
    nifti_image * nii_columnarity = nifti_copy_nim_info(nii);
    nii_columnarity->nt = 1;
    nii_columnarity->nvox = nii->nvox / size_time;
    nii_columnarity->datatype = NIFTI_TYPE_FLOAT32;
    nii_columnarity->nbyper = sizeof(float);
    nii_columnarity->data = calloc(nii_columnarity->nvox, nii_columnarity->nbyper);
    float* nii_columnarity_data = static_cast<float*>(nii_columnarity->data);

    nifti_image * nii_parcelval = nifti_copy_nim_info(nii);
    nii_parcelval->nt = 1;
    nii_parcelval->nvox = nii->nvox / size_time;
    nii_parcelval->datatype = NIFTI_TYPE_FLOAT32;
    nii_parcelval->nbyper = sizeof(float);
    nii_parcelval->data = calloc(nii_parcelval->nvox, nii_parcelval->nbyper);
    float* nii_parcelval_data = static_cast<float*>(nii_parcelval->data);

    // ========================================================================
    ////////////////////////////////////////////////////////////
    // finding number of layers and columsn to work with ///////
    ////////////////////////////////////////////////////////////
    
    int Nr_columns = 0; 
    int Nr_layers  = 0; 
    int Nr_parcels  = 0; 
    
    for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		   if ((int)*(nim_columns_data  +  nxy*iz + nx*iy + ix) >= Nr_columns ) {
			    Nr_columns = (int)*(nim_columns_data +  nxy*iz + nx*iy + ix) ; 
		   }
		   if ((int)*(nim_layers_data   +  nxy*iz + nx*iy + ix) >= Nr_layers ) {
			    Nr_layers = (int)*(nim_layers_data +  nxy*iz + nx*iy + ix) ; 
		   }
	    }
	  }
    }
    cout << "there are "<< Nr_columns << " columns"   << endl;
    cout << "there are "<< Nr_layers  << "  layers" << endl << endl;
    Nr_parcels = Nr_columns  * Nr_layers ; 
    cout << "there are "<< Nr_parcels  << " parcels" << endl << endl;
    
    
    // ========================================================================
    /////////////////////////////////////////////////////////////////////
    // finding approximate number of voxels per parcel (layer column combination) ///////
    /////////////////////////////////////////////////////////////////////
    int layeridx = 0; 
    int columnindx = 0; 
    int voxelidx = 0;
    int runningidx = 0; 
    
    double vec_nrVox_pacels[Nr_parcels] ; // Access the parcel as; vec_nrVox_pacels[ Nr_columns * LayerIndex + ColumnIndex]
    double vec_mostcommonval_pacels[Nr_parcels] ; // Access the parcel as; vec_mostcommonval_pacels[ Nr_columns * LayerIndex + ColumnIndex]

    
    for (int ip=0; ip<Nr_parcels; ++ip){
			vec_nrVox_pacels[ip] = 0.; 
			vec_mostcommonval_pacels[ip] = 0.; 
	}
	
	for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		 //  vec_val_pacels[(int)*(nim_columns_data +  nxy*iz + nx*iy + ix)][(int)*(nim_layers_data +  nxy*iz + nx*iy + ix)] =+  ;
		 //cout <<  Nr_columns * (int)*(nim_layers_data +  nxy*iz + nx*iy + ix) + (int)*(nim_columns_data +  nxy*iz + nx*iy + ix) << endl; 
		 layeridx = (int)*(nim_layers_data +  nxy*iz + nx*iy + ix)-1;
		 columnindx = (int)*(nim_columns_data +  nxy*iz + nx*iy + ix)-1; 
		 vec_nrVox_pacels[ Nr_columns * layeridx + columnindx]++ ; 
	    }
	  }
    }
        
    double mean_Nr_voxels_per_parcel = 0; 
    double stdev_Nr_voxels_per_parcel = 0; 
    
    mean_Nr_voxels_per_parcel =  ren_average(vec_nrVox_pacels, Nr_layers*Nr_columns); 
    stdev_Nr_voxels_per_parcel = ren_stdev(vec_nrVox_pacels, Nr_layers*Nr_columns); 
    
    cout << "Mean number of voxels per pacel is " << mean_Nr_voxels_per_parcel << endl; 
    cout << "Stdev number of voxels per pacel is " << stdev_Nr_voxels_per_parcel << endl;   
    
    
    // ========================================================================
    /////////////////////////////////////////////////////////////////////
    // finding most common value in each parcel ///////
    /////////////////////////////////////////////////////////////////////   
    
    int MaxVoxelPerParcel = (int) mean_Nr_voxels_per_parcel + (int) stdev_Nr_voxels_per_parcel; // for very very large columns, some voxels are missed.
      cout << "MaxVoxelPerParcel is " << MaxVoxelPerParcel << endl;   

    int* val_Vox_pacels = new int[Nr_parcels*MaxVoxelPerParcel]; 
    //Access the parcel as: val_Vox_pacels[ MaxVoxelPerParcel * Nr_columns * LayerIndex + MaxVoxelPerParcel * ColumnIndex + VoxelIndex]
    
    
    /////////////////////////////////////////////////////////////////////
    // loop across voxels for filling vector across parcels  ///////
    /////////////////////////////////////////////////////////////////////   
    
    for (int il=0; il<Nr_layers; ++il){
		for (int ic=0; ic<Nr_columns; ++ic){	
			    vec_nrVox_pacels[Nr_columns * il +  ic] = 0.; 
			    
			   for (int iv=0; iv<MaxVoxelPerParcel; ++iv){
			    val_Vox_pacels[ MaxVoxelPerParcel * Nr_columns * il + MaxVoxelPerParcel * ic + iv] = 0 ; 
		       }
		}
	}
                  cout << " Hi1 "  << endl;   

    
   for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		   if (*(nim_columns_data +  nxy*iz + nx*iy + ix) >0) {     // NOTE(Renzo): include a check somewhere above that the layer-file makes sense in all voxles where columns are defined.
			layeridx = *(nim_layers_data +  nxy*iz + nx*iy + ix) -1 ;
			columnindx = *(nim_columns_data +  nxy*iz + nx*iy + ix) -1 ; 
			vec_nrVox_pacels[Nr_columns * layeridx +  columnindx] += 1 ; 
			voxelidx = vec_nrVox_pacels[Nr_columns * layeridx +  columnindx] - 1 ; 
			voxelidx = min(voxelidx,MaxVoxelPerParcel-1);  // ignore the last voxels of gigantic column outliers. 
		 
			//cout << layeridx << " " << columnindx << "  "<< voxelidx << endl ; 
			val_Vox_pacels[MaxVoxelPerParcel * Nr_columns * layeridx + MaxVoxelPerParcel * columnindx + voxelidx] =  *(nii_data +  nxy*iz + nx*iy + ix) ; 
		   }
	   }
	 }
   }
   
    ///////////////////////////////////////////////////////////
    // finding most common value in each parcel output  ///////
    ///////////////////////////////////////////////////////////
    
    cout << " Hi1a "  << endl;  
    //vec_mostcommonval_pacels[ Nr_columns * LayerIndex + ColumnIndex]
    
    //double vec_val_in_parcel[MaxVoxelPerParcel] ; 
    int* vec_val_in_parcel = new int[MaxVoxelPerParcel];
    for (int ivox=0; ivox<MaxVoxelPerParcel; ++ivox) vec_val_in_parcel[ivox] = 0; 
	int NRvoxel_in_this_parcel = 0; 
   

    for (int il=0; il<Nr_layers; ++il){
		for (int ic=0; ic<Nr_columns; ++ic){
			  
			    NRvoxel_in_this_parcel = vec_nrVox_pacels[Nr_columns * il +  ic]; 
			    NRvoxel_in_this_parcel = min(MaxVoxelPerParcel,NRvoxel_in_this_parcel); 
			    
			  
				for (int ivox=0; ivox<NRvoxel_in_this_parcel; ++ivox){
                    vec_val_in_parcel[ivox] =  val_Vox_pacels[MaxVoxelPerParcel * Nr_columns * il + MaxVoxelPerParcel * ic + ivox]; 
			    }
			    
			    vec_mostcommonval_pacels[ Nr_columns * il + ic] = ren_most_occurred_number(vec_val_in_parcel, NRvoxel_in_this_parcel); 
			    
			   // cout << il << " " << ic << "  "<< vec_mostcommonval_pacels[ Nr_columns * il + ic]<< endl;
		}
	}
    

   for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		   if (*(nim_columns_data +  nxy*iz + nx*iy + ix) >0) { 
			layeridx = *(nim_layers_data +  nxy*iz + nx*iy + ix) -1 ;
			columnindx = *(nim_columns_data +  nxy*iz + nx*iy + ix) -1 ; 
		    *(nii_parcelval_data +  nxy*iz + nx*iy + ix) = vec_mostcommonval_pacels[ Nr_columns * layeridx + columnindx];
		   }
	   }
	 }
   }
   
    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "parcel_val", nii_parcelval, true, use_outpath);
   
    ////////////////////////////////////////////////////////////////////////////////////
    // finding local patch of neigbors (snipets repurposed from LN2_COLTRAIN) ///////
    ////////////////////////////////////////////////////////////////////////////////////
    uint32_t ix, iy, iz;
    uint32_t iix, iiy, iiz;
	// Voxels of interest
    uint32_t nr_voi = 0;
    for (uint32_t i = 0; i != nxyz; ++i) {
        if (*(nim_columns_data + i) > 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nxyz; ++i) {
        if (*(nim_columns_data + i) > 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Triangulation on subset Voronoi cells
    // ========================================================================
    // Keep neighbouring column ID's in an array
   int32_t* neighbours;
   neighbours = (int32_t*) malloc(27*sizeof(int32_t));// there are 26 neighbors of in a rubics cube.
   for (int j=0; j<27; ++j) {
        *(neighbours + j) = 0;
   }
    
   int* ColId_in_neigh_column = new int[Nr_parcels*27];
      //Access the parcel as: [ Nr_parcels * NeighColidx  +   parcelID ]  
   // center parcelID: Nr_columns * layeridx + columnindx
   int* val_in_neigh_column   = new int[Nr_parcels*27];
      //Access the parcel as: [ Nr_parcels * NeighColidx +   parcelID ]  
   // center parcelID: Nr_columns * layeridx + columnindx
   for (int j=0; j<Nr_parcels*27; ++j) {
	   ColId_in_neigh_column [j] =0 ;
	   val_in_neigh_column [j] =0 ;
   }
   
   int* Collumn_in_neighbor = new int[Nr_parcels*27];
   int* layer_in_neighbor   = new int[Nr_parcels*27];
   for (int j=0; j<Nr_parcels*27; ++j) {
	   Collumn_in_neighbor [j] =0 ;
	   layer_in_neighbor [j] =0 ;
   }
      

   int runidx_neigh_parcel[Nr_parcels];
      //Access the parcel as: [ parcelID ]  
   // center parcelID: Nr_columns * layeridx + columnindx
   for (int j=0; j<Nr_parcels; ++j) {
	   runidx_neigh_parcel [j] = 1 ;
   }
   
  
   
   int current_column_idx = 0 ;
   int current_layer_idx  = 0 ; 
   int current_parcel_idx  = 0 ; 
   int current_column_neighbor = 0 ;
   int current_layer_neighbor  = 0 ; 
   int current_parcel_neighbor  = 0 ;
   int layer_of_new_neigbor = 0; 
   int niegbor_id = 0; 
   int is_notnew_nighbor = 0;
   int is_notnew_parcel_nighbor = 0;
//   int is_notnew_layer_nighbor = 0;
   
    cout << " Hi3 "  << endl;  
  
    for(uint32_t ii = 0; ii != nr_voi; ++ii) { 
        uint32_t i = *(voi_id + ii);
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        current_column_idx = (int)*(nim_columns_data + i) -1;
        current_layer_idx  = (int)*(nim_layers_data  + i) -1;
        current_parcel_idx =  Nr_columns * current_layer_idx + current_column_idx;
        
        //cout << "column, layer, parcel " << current_column_idx << " "<< current_layer_idx << " "<< current_parcel_idx << " "<< endl;
        //*(neighbours + 0) = *(nim_columns_data + i);

		layer_in_neighbor   	[ Nr_parcels * 0 +  current_parcel_idx ] = current_layer_idx ; 
        Collumn_in_neighbor 	[ Nr_parcels * 0 +  current_parcel_idx ] = current_column_idx ; 
        val_in_neigh_column 	[ Nr_parcels * 0 +  current_parcel_idx ] = (int)*(nii_parcelval_data + i) ;
        ColId_in_neigh_column 	[ Nr_parcels * 0 +  current_parcel_idx ] = current_parcel_idx; 
        
        // 1-jump neighbours
         if (ix > 0) {
            int j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
				
				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
				 
					 
//				if (current_parcel_neighbor !=current_parcel_idx &&   is_notnew_parcel_nighbor == 0 ) {
//					 cout << " current_parcel_idx " << current_parcel_idx << " current_column_idx " << current_column_idx << " current_layer_idx " << current_layer_idx  << endl;
//					 cout << " current_parcel_neighbor " << current_parcel_neighbor << " current_column_neighbor " << current_column_neighbor << " current_layer_neighbor " << current_layer_neighbor  << endl;
//				     cout << ix<< " "<< iy << " " << iz <<"       " <<   runidx_neigh_parcel  [current_parcel_idx]<< "  "<< runidx_neigh_parcel  [current_parcel_idx]<<  endl ; 
//				     cout << ix-1<< " "<< iy << " " << iz <<"       " <<   runidx_neigh_parcel  [current_parcel_idx]<< "  "<< runidx_neigh_parcel  [current_parcel_idx]<<  endl ; 

				
//					      for (int niegbor_id_index = 0 ; niegbor_id_index < 27; niegbor_id_index++ ){ 
//					      	cout << ColId_in_neigh_column[Nr_parcels * niegbor_id_index +  current_parcel_idx ] << " " ;
					      	//cout << ColId_in_neigh_layer [Nr_parcels * current_column_neighboridx +  current_parcel_idx ] << " " ;
//						  }
//						  cout << endl; 
//				          for (int niegbor_id_index = 0 ; niegbor_id_index < 27; niegbor_id_index++ ){ 
//					      	cout << Collumn_in_neighbor[Nr_parcels * niegbor_id_index +  current_parcel_idx ] << " " ;
//				          }     
//				          cout << endl; 
//				          for (int niegbor_id_index = 0 ; niegbor_id_index < 27; niegbor_id_index++ ){ 
//					      	cout << layer_in_neighbor[Nr_parcels * niegbor_id_index +  current_parcel_idx ] << " " ;
//				          }     
//				          cout << endl; 
//				          for (int niegbor_id_index = 0 ; niegbor_id_index < 27; niegbor_id_index++ ){ 
//					      	cout << val_in_neigh_column[Nr_parcels * niegbor_id_index +  current_parcel_idx ] << " " ;
					      	//cout << ColId_in_neigh_layer [Nr_parcels * current_column_neighboridx +  current_parcel_idx ] << " " ;
//						  }
//						  cout << endl; 
//				}
				    
            } //loop across neigbors closed
        }
        
        		
				 
        if (iy > 0) {
            int j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
               				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iz > 0) {
            int j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix < size_x) {
            int j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iy < size_y) {
            int j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iz < size_z) {
            int j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }

        // 2-jump neighbours
        if (ix > 0 && iy > 0) {
            int j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix > 0 && iy < size_y) {
            int j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix < size_x && iy > 0) {
            int j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix < size_x && iy < size_y) {
            int j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
       if (iy > 0 && iz > 0) {
            int j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iy > 0 && iz < size_z) {
            int j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iy < size_y && iz > 0) {
            int j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (iy < size_y && iz < size_z) {
            int j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix > 0 && iz > 0) {
            int j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
               				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix < size_x && iz > 0) {
            int j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
               				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix > 0 && iz < size_z) {
            int j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }
        if (ix < size_x && iz < size_z) {
            int j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                				 is_notnew_parcel_nighbor = 0; 
				 current_column_neighbor  = (int)*(nim_columns_data + j)-1;
                 current_layer_neighbor   = (int)*(nim_layers_data + j) -1; 
                 current_parcel_neighbor  =  Nr_columns * current_layer_neighbor + current_column_neighbor ;
                 
               // checking if this neigbor is already considered or not        
				 for (int current_column_neighboridx = 0 ; current_column_neighboridx < 27; current_column_neighboridx++ ){ 
					     if (current_parcel_neighbor == ColId_in_neigh_column[Nr_parcels * current_column_neighboridx +  current_parcel_idx ]){
					         is_notnew_parcel_nighbor ++;
						 }
				 }
				 if (current_parcel_idx == current_parcel_neighbor )   is_notnew_parcel_nighbor = 1;
			 
				 if ( is_notnew_parcel_nighbor == 0 ){
					 niegbor_id = min (runidx_neigh_parcel [current_parcel_idx], 27) ; // this is to make sure that we are not runnign out of memory if there are more than 27 nioghbors. 
					 ColId_in_neigh_column [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_parcel_neighbor;
					 val_in_neigh_column   [Nr_parcels * niegbor_id +  current_parcel_idx ] = (int)*(nii_parcelval_data + j) ; 
					 Collumn_in_neighbor   [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_column_neighbor ; 
					 layer_in_neighbor     [Nr_parcels * niegbor_id +  current_parcel_idx ] = current_layer_neighbor ; 
					 runidx_neigh_parcel   [current_parcel_idx] ++; 
				 }	 
            }
        }

        //  not useing 3-jump neighbours 
        if (ix > 0 && iy > 0 && iz > 0) {
            int j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 19) = *(nim_columns_data + j);
            }
        }
        if (ix > 0 && iy > 0 && iz < size_z) {
            int j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
               // *(neighbours + 20) = *(nim_columns_data + j);
            }
        }
        if (ix > 0 && iy < size_y && iz > 0) {
            int j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 21) = *(nim_columns_data + j);
            }
        }
        if (ix < size_x && iy > 0 && iz > 0) {
            int j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 22) = *(nim_columns_data + j);
            }
        }
        if (ix > 0 && iy < size_y && iz < size_z) {
            int j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 23) = *(nim_columns_data + j);
            }
        }
        if (ix < size_x && iy > 0 && iz < size_z) {
            int j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 24) = *(nim_columns_data + j);
            }
        }
        if (ix < size_x && iy < size_y && iz > 0) {
            int j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 25) = *(nim_columns_data + j);
            }
        }
        if (ix < size_x && iy < size_y && iz < size_z) {
            int j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
            if (*(nim_columns_data + j) > 0) {
                //*(neighbours + 26) = *(nim_columns_data + j);
            }
        }



		/*
        // Have an array to note down columns id's that will form triangles.
        // NOTE(Faruk): I allocate 8 elements here because if the volume is
        // over-parcellated (i.e. too many columns), hit count can go above 3
        int32_t* trio;
        trio = (int32_t*) malloc(8*sizeof(int32_t));

        for (int m=0; m<8; m++) {
            int n;
            for (n=0; n<m; n++) {
                if (*(neighbours + m) == *(neighbours + n)) {
                    break;
                }
            }
            if (m == n) {
                if (*(neighbours + m) != 0) {
                    *(trio + count_neigh) = *(neighbours + m);
                    count_neigh += 1;
                }
            }
        }
        if (count_neigh == 3) {  // Hard-lock to triplets only
            for (int l=0; l<count_neigh; ++l) {
                *(triplets + count_triplets * 3 + l) = *(trio + l);
            }
            count_triplets += 1;
        }
        */
    }



cout << " I am done filling the neighbor array "<< endl; 

//////////////////////////////////////////////////////////////
////////////  Loop over parcels and quantifying their similary of valued across layers and columns
//////////////////////////////////////////////////////////////////////

double laminarity[Nr_parcels];      //Access the parcel as: [ parcelID ]  // center parcelID: Nr_columns * layeridx + columnindx
for (int j=0; j<Nr_parcels; ++j)   laminarity [j] = 0;  
int count_same_val_layer = 0;    
int count_same_layer = 0;       

for (int paridx=0; paridx<Nr_parcels; ++paridx) {
//	cout << "runidx_neigh_parcel["<<paridx<<"]   "<< runidx_neigh_parcel[paridx] << endl; 
//    cout << "   layer_in_neighbor[Nr_parcels *"<<paridx<<" + 0]   "<< layer_in_neighbor[Nr_parcels * paridx +  0 ] << endl; 
	count_same_val_layer = 0;    
	count_same_layer = 0;   

	for (int nigidx=1; nigidx<min(runidx_neigh_parcel[paridx],27); ++nigidx){
		if (layer_in_neighbor[Nr_parcels * 0 +  paridx] == layer_in_neighbor[Nr_parcels * nigidx +  paridx ]  ){
			count_same_layer ++;
			if (val_in_neigh_column[Nr_parcels * 0 +  paridx ] == val_in_neigh_column[Nr_parcels * nigidx +  paridx ]  ){
				count_same_val_layer++;
			}
		}
	}	
	if ( count_same_layer > 0 ) laminarity [paridx] = (double)count_same_val_layer / (double) count_same_layer; 
	 //cout << "   laminarity[ "<<paridx<<" ]   "<< laminarity [paridx]<< endl; 
}

   for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		   if (*(nim_columns_data +  nxy*iz + nx*iy + ix) >0) { 
			layeridx = *(nim_layers_data +  nxy*iz + nx*iy + ix) -1 ;
			columnindx = *(nim_columns_data +  nxy*iz + nx*iy + ix) -1 ; 
		    *(nii_laminarity_data +  nxy*iz + nx*iy + ix)  = laminarity [ Nr_columns * layeridx + columnindx];
		   }
	   }
	 }
   }


    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "output_laminarity", nii_laminarity, true, use_outpath);
    

double columnarity[Nr_parcels];      //Access the parcel as: [ parcelID ]  // center parcelID: Nr_columns * layeridx + columnindx
for (int j=0; j<Nr_parcels; ++j)   laminarity [j] = 0;  
int count_same_val_col = 0;    
int count_same_col = 0;       

for (int paridx=0; paridx<Nr_parcels; ++paridx) {
//	cout << "runidx_neigh_parcel["<<paridx<<"]   "<< runidx_neigh_parcel[paridx] << endl; 
//    cout << "   layer_in_neighbor[Nr_parcels *"<<paridx<<" + 0]   "<< layer_in_neighbor[Nr_parcels * paridx +  0 ] << endl; 
	count_same_val_col = 0;    
	count_same_col = 0;   

	for (int nigidx=1; nigidx<min(runidx_neigh_parcel[paridx],27); ++nigidx){
		if (Collumn_in_neighbor[Nr_parcels * 0 +  paridx] == Collumn_in_neighbor[Nr_parcels * nigidx +  paridx ]  ){
			count_same_col ++;
			if (val_in_neigh_column[Nr_parcels * 0 +  paridx ] == val_in_neigh_column[Nr_parcels * nigidx +  paridx ]  ){
				count_same_val_col++;
			}
		}
	}	
	if ( count_same_col > 0 ) columnarity [paridx] = (double)count_same_val_col / (double) count_same_col; 
	 //cout << "   columnarity[ "<<paridx<<" ]   "<< columnarity [paridx]<< endl; 
}


   for(int iy=0; iy<size_y; ++iy){
     for(int ix=0; ix<size_x; ++ix){
       for(int iz=0; iz<size_z; ++iz){
		   if (*(nim_columns_data +  nxy*iz + nx*iy + ix) >0) { 
			layeridx = *(nim_layers_data +  nxy*iz + nx*iy + ix) -1 ;
			columnindx = *(nim_columns_data +  nxy*iz + nx*iy + ix) -1 ; 
		    *(nii_laminarity_data +  nxy*iz + nx*iy + ix)  = laminarity [ Nr_columns * layeridx + columnindx];
		    *(nii_columnarity_data +  nxy*iz + nx*iy + ix) = columnarity[ Nr_columns * layeridx + columnindx];
		   }
	   }
	 }
   }

//Collumn_in_neighbor[Nr_parcels * niegbor_id_index +  current_parcel_idx ] 
//layer_in_neighbor[Nr_parcels * niegbor_id_index +  current_parcel_idx ] 
//val_in_neigh_column[Nr_parcels * niegbor_id_index +  current_parcel_idx ] 
//runidx_neigh_parcel   [current_parcel_idx]



    // =======================
    //////////////////////////
    // writing output  ///////
    //////////////////////////

    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "output_columnarity", nii_columnarity, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}

  
