# For internal testing. Just to check whether programs execute.

./LN_LAYER_SMOOTH -input UNI.nii -layer_file layers_rim.nii -FWHM 1
./LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 40 -shift
./LN_MP2RAGE_DNOISE -INV1 INV1.nii -INV2 INV2.nii -UNI UNI.nii

# ./LN_3DCOLUMNS -layer_file layers_large.nii -landmarks landmarks.nii
./LN_CORREL2FILES -file1 lo_Nulled_intemp.nii -file2 lo_BOLD_intemp.nii
# ./LN_COLUMNAR_DIST
./LN_DIRECT_SMOOTH -input UNI.nii -FWHM 2 -direction 1
./LN_EXTREMETR -input lo_BOLD_intemp.nii
./LN_FLOAT_ME -input rim.nii
./LN_GFACTOR -input activity_map_example.nii -output gfactor -variance 1 -direction 1 -grappa 2 -cutoff 5
./LN_GRADSMOOTH -gradfile VASO_LN_Tmin.nii -input VASO_LN_Tmean.nii -FWHM 2 -within -selectivity 0.1
./LN_GROW_LAYERS -rim rim.nii
./LN_INTPRO -image UNI.nii -min -direction 2

# ./LN_LEAKY_LAYERS -rim rim_test.nii
./LN_NOISEME -input activity_map_example.nii -output noiseme_custom_output
# ./LN_RAGRUG -input activity_map_example.nii
./LN_SKEW -timeseries lo_BOLD_intemp.nii
./LN_TEMPSMOOTH -input VASO_LN.nii -box 1
./LN_TEMPSMOOTH -input VASO_LN.nii -gaus 1
./LN_TRIAL -file BOLD_intemp.nii -trialdur 12
./LN_ZOOM -mask UNI_MASK.nii -input UNI.nii
