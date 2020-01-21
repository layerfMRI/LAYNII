#!/bin/bash

g++    -c -o nifti2_io.o nifti2_io.cpp
g++    -c -o nifticdf.o nifticdf.cpp
g++    -c -o znzlib.o znzlib.cpp
g++    -c -o LN_FAsim.o LN_FAsim.cpp
g++  -o LN_FAsim -Wall -pedantic -DHAVE_ZLIB -I.  LN_FAsim.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_NOISEME.o LN_NOISEME.cpp
g++  -o LN_NOISEME -Wall -pedantic -DHAVE_ZLIB -I.  LN_NOISEME.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_GROW_LAYERS.o LN_GROW_LAYERS.cpp
g++  -o LN_GROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  LN_GROW_LAYERS.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_3DGROW_LAYERS.o LN_3DGROW_LAYERS.cpp
g++  -o LN_3DGROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  LN_3DGROW_LAYERS.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_DEBUGGING.o LN_DEBUGGING.cpp
g++  -o LN_DEBUGGING -Wall -pedantic -DHAVE_ZLIB -I.  LN_DEBUGGING.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_GFACTOR.o LN_GFACTOR.cpp
g++  -o LN_GFACTOR -Wall -pedantic -DHAVE_ZLIB -I.  LN_GFACTOR.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_LEAKY_LAYERS.o LN_LEAKY_LAYERS.cpp
g++  -o LN_LEAKY_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  LN_LEAKY_LAYERS.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_LAYER_SMOOTH.o LN_LAYER_SMOOTH.cpp
g++  -o LN_LAYER_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  LN_LAYER_SMOOTH.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_3DCOLUMNS.o LN_3DCOLUMNS.cpp
g++  -o LN_3DCOLUMNS -Wall -pedantic -DHAVE_ZLIB -I.  LN_3DCOLUMNS.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_SHORT_ME.o LN_SHORT_ME.cpp
g++  -o LN_SHORT_ME -Wall -pedantic -DHAVE_ZLIB -I.  LN_SHORT_ME.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_FIX_RIM.o LN_FIX_RIM.cpp
g++  -o LN_FIX_RIM -Wall -pedantic -DHAVE_ZLIB -I.  LN_FIX_RIM.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_FLOAT_ME.o LN_FLOAT_ME.cpp
g++  -o LN_FLOAT_ME -Wall -pedantic -DHAVE_ZLIB -I.  LN_FLOAT_ME.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_IMAGIRO.o LN_IMAGIRO.cpp
g++  -o LN_IMAGIRO -Wall -pedantic -DHAVE_ZLIB -I.  LN_IMAGIRO.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_DIRECT_SMOOTH.o LN_DIRECT_SMOOTH.cpp
g++  -o LN_DIRECT_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  LN_DIRECT_SMOOTH.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_RAGRUG.o LN_RAGRUG.cpp
g++  -o LN_RAGRUG -Wall -pedantic -DHAVE_ZLIB -I.  LN_RAGRUG.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_CORREL2FILES.o LN_CORREL2FILES.cpp
g++  -o LN_CORREL2FILES -Wall -pedantic -DHAVE_ZLIB -I.  LN_CORREL2FILES.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_EXTREMETR.o LN_EXTREMETR.cpp
g++  -o LN_EXTREMETR -Wall -pedantic -DHAVE_ZLIB -I.  LN_EXTREMETR.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_BOCO.o LN_BOCO.cpp
g++  -o LN_BOCO -Wall -pedantic -DHAVE_ZLIB -I.  LN_BOCO.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_TRIAL.o LN_TRIAL.cpp
g++  -o LN_TRIAL -Wall -pedantic -DHAVE_ZLIB -I.  LN_TRIAL.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_ZOOM.o LN_ZOOM.cpp
g++  -o LN_ZOOM -Wall -pedantic -DHAVE_ZLIB -I.  LN_ZOOM.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_SMOOTH_RIM.o LN_SMOOTH_RIM.cpp
g++  -o LN_SMOOTH_RIM -Wall -pedantic -DHAVE_ZLIB -I.  LN_SMOOTH_RIM.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_COLUMNAR_DIST.o LN_COLUMNAR_DIST.cpp
g++  -o LN_COLUMNAR_DIST -Wall -pedantic -DHAVE_ZLIB -I.  LN_COLUMNAR_DIST.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_GRADSMOOTH.o LN_GRADSMOOTH.cpp
g++  -o LN_GRADSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  LN_GRADSMOOTH.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_SKEW.o LN_SKEW.cpp
g++  -o LN_SKEW -Wall -pedantic -DHAVE_ZLIB -I.  LN_SKEW.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_INTPRO.o LN_INTPRO.cpp
g++  -o LN_INTPRO -Wall -pedantic -DHAVE_ZLIB -I.  LN_INTPRO.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_TEMPSMOOTH.o LN_TEMPSMOOTH.cpp
g++  -o LN_TEMPSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  LN_TEMPSMOOTH.o nifti2_io.o nifticdf.o znzlib.o
g++    -c -o LN_MP2RAGE_DNOISE.o LN_MP2RAGE_DNOISE.cpp
g++  -o LN_MP2RAGE_DNOISE -Wall -pedantic -DHAVE_ZLIB -I.  LN_MP2RAGE_DNOISE.o nifti2_io.o nifticdf.o znzlib.o 
g++    -c -o LN_PHYSIO_PARS.o LN_PHYSIO_PARS.cpp
g++  -o LN_PHYSIO_PARS -Wall -pedantic -DHAVE_ZLIB -I.  LN_PHYSIO_PARS.o nifti2_io.o nifticdf.o znzlib.o 
