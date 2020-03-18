# For internal testing. Just to check whether programs execute.

../LN_LAYER_SMOOTH -input sc_VASO_stat.nii -layer_file sc_layers.nii -FWHM 1
../LN_BOCO -Nulled lo_Nulled_intemp.nii -BOLD lo_BOLD_intemp.nii -trialBOCO 40 -shift
../LN_MP2RAGE_DNOISE -INV1 sc_INV1.nii -INV2 sc_INV2.nii -UNI sc_UNI.nii

../LN_3DCOLUMNS -layers sc_layers_3dcolumns.nii -landmarks sc_landmarks_3dcolumns.nii
../LN_CORREL2FILES -file1 lo_Nulled_intemp.nii -file2 lo_BOLD_intemp.nii
../LN_COLUMNAR_DIST -layers sc_layers_3dcolumns.nii -landmarks sc_landmarks_3dcolumns.nii
../LN_DIRECT_SMOOTH -input sc_UNI.nii -FWHM 2 -direction 1
../LN_EXTREMETR -input lo_BOLD_intemp.nii
../LN_FLOAT_ME -input sc_rim.nii
../LN_GFACTOR -input sc_UNI.nii -output gfactor -variance 1 -direction 1 -grappa 2 -cutoff 5
../LN_GRADSMOOTH -gradfile lo_T1EPI.nii -input lo_VASO_stat.nii -FWHM 2 -within -selectivity 0.1
../LN_GROW_LAYERS -rim sc_rim.nii
../LN_INTPRO -image sc_UNI.nii -min -direction 2 -range 3

../LN_LEAKY_LAYERS -rim sc_rim.nii
../LN_NOISEME -input lo_VASO_stat.nii -std 1
../LN_RAGRUG -input sc_rim.nii
../LN_SKEW -input lo_BOLD_intemp.nii
../LN_TEMPSMOOTH -input lo_BOLD_intemp.nii -box 1
../LN_TEMPSMOOTH -input lo_BOLD_intemp.nii -gaus 1
../LN_TRIAL -input lo_BOLD_intemp.nii -trialdur 20
../LN_ZOOM -mask sc_layers_3dcolumns.nii -input sc_UNI.nii
