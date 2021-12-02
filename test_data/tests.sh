# For internal testing. Just to check whether programs execute.

../LN2_LAYER_SMOOTH -input sc_VASO_act.nii.gz -layer_file sc_layers.nii.gz -FWHM 1
../LN_LAYER_SMOOTH -input sc_VASO_act.nii.gz -layer_file sc_layers.nii.gz -FWHM 0.3 -NoKissing

../LN_BOCO -Nulled lo_Nulled_intemp.nii.gz -BOLD lo_BOLD_intemp.nii.gz -trialBOCO 40 -shift
../LN_MP2RAGE_DNOISE -INV1 sc_INV1.nii.gz -INV2 sc_INV2.nii.gz -UNI sc_UNI.nii.gz

../LN2_LAYERS -rim sc_rim.nii.gz -nr_layers 10 -equivol

../LN_3DCOLUMNS -layers sc_layers_3dcolumns.nii.gz -landmarks sc_landmarks_3dcolumns.nii.gz
../LN_CORREL2FILES -file1 lo_Nulled_intemp.nii.gz -file2 lo_BOLD_intemp.nii.gz
../LN_COLUMNAR_DIST -layers sc_layers_3dcolumns.nii.gz -landmarks sc_landmarks.nii.gz
../LN_DIRECT_SMOOTH -input sc_UNI.nii.gz -FWHM 2 -direction 3
../LN_EXTREMETR -input lo_BOLD_intemp.nii.gz
../LN_FLOAT_ME -input lo_BOLD_intemp.nii.gz
../LN_SHORT_ME -input lo_VASO_act.nii.gz -output short.nii.gz
../LN_INT_ME -input LN_INT_ME -input lo_BOLD_act.nii.gz
../LN_GFACTOR -input sc_INV2.nii.gz  -variance 1 -direction 1 -grappa 2 -cutoff 200
../LN_GRADSMOOTH -gradfile lo_gradT1.nii.gz -input lo_VASO_act.nii.gz -FWHM 1 -within -selectivity 0.1
../LN_GRADSMOOTH_ITER -gradfile lo_gradT1.nii.gz -input lo_VASO_act.nii.gz -FWHM 1 -within -selectivity 0.1
../LN_GROW_LAYERS -rim sc_rim.nii.gz
../LN_INTPRO -image sc_UNI.nii.gz -min -direction 2 -range 3
../LN_IMAGIRO -layers sc_layers_3dcolumns.nii.gz -columns sc_columns_3dcolumns.nii.gz -data sc_BOLD_act.nii.gz 
../LN_LEAKY_LAYERS -rim lo_rim_LL.nii.gz
../LN_NOISEME -input lo_VASO_act.nii.gz -std 1
../LN_RAGRUG -input sc_rim.nii.gz
../LN_SKEW -input lo_BOLD_intemp.nii.gz
../LN_TEMPSMOOTH -input lo_BOLD_intemp.nii.gz -box 1
../LN_TEMPSMOOTH -input lo_BOLD_intemp.nii.gz -gaus 1
../LN_TRIAL -input lo_BOLD_intemp.nii.gz -trialdur 20
../LN_ZOOM -mask sc_layers_3dcolumns.nii.gz -input sc_UNI.nii.gz
../LN_LOITUMA -equidist sc_distlay_1000.nii.gz -leaky sc_leakylay_1000.nii.gz -FWHM 1 -nr_layers 10
../LN_NOISE_KERNEL -input lo_Nulled_intemp.nii.gz -kernel_size 7
../LN2_DEVEIN -layer_file lo_layers.nii.gz -column_file lo_columns.nii.gz -input lo_BOLD_act.nii.gz -ALF lo_ALF.nii.gz
../LN2_RIMIFY -input sc_rim.nii.gz -innergm 2 -outergm 1 -gm 3 -output rimified_tim.nii.gz
../LN_INFO -input lo_T1EPI.nii.gz
../LN_CONLAY -layers lo_sc_layers.nii.gz -ref lo_T1EPI.nii.gz -subsample -output lo_layers_out.nii.gz
../LN2_COLUMNS -rim sc_rim.nii.gz -midgm sc_midGM.nii.gz -nr_columns 300
../LN2_CHOLMO -layers sc_layers.nii.gz -outer -nr_layers 3 -layer_thickness 0.4 -output padded_layers.nii.gz
../LN2_PROFILE -input sc_VASO_act.nii.gz -layers sc_layers.nii.gz -plot
../LN2_LAYERDIMENSION -values lo_BOLD_act.nii.gz -layers lo_layers.nii.gz -columns lo_columns.nii.gz
../LN2_MASK -scores lo_BOLD_act.nii.gz -columns lo_columns.nii.gz -mean_thr 1 -output mask.nii.gz -abs
