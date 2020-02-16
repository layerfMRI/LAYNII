#!/bin/bash

# Dependencies
g++    -c -o obj/nifti2_io.o deps/nifti2_io.cpp
g++    -c -o obj/nifticdf.o deps/nifticdf.cpp
g++    -c -o obj/znzlib.o deps/znzlib.cpp
g++    -c -o obj/renzo_stat.o src/renzo_stat.c
g++    -c -o obj/utils.o src/utils.cpp

# # LAYNII
g++    -c -o obj/LN_MP2RAGE_DNOISE.o src/LN_MP2RAGE_DNOISE.cpp
g++  -o bin/LN_MP2RAGE_DNOISE -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_MP2RAGE_DNOISE.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/utils.o
g++    -c -o obj/LN_BOCO.o src/LN_BOCO.cpp
g++  -o bin/LN_BOCO -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_BOCO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/renzo_stat.o obj/utils.o
g++    -c -o obj/LN_LAYER_SMOOTH.o src/LN_LAYER_SMOOTH.cpp
g++  -o bin/LN_LAYER_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_LAYER_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/utils.o
# g++    -c -o obj/LN_CORREL2FILES.o src/LN_CORREL2FILES.cpp
# g++  -o bin/LN_CORREL2FILES -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_CORREL2FILES.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/renzo_stat.o
# g++    -c -o obj/LN_GRADSMOOTH.o src/LN_GRADSMOOTH.cpp
# g++  -o bin/LN_GRADSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GRADSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/renzo_stat.o
# g++    -c -o obj/LN_SKEW.o src/LN_SKEW.cpp
# g++  -o bin/LN_SKEW -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SKEW.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/renzo_stat.o
# g++    -c -o obj/LN_FAsim.o src/LN_FAsim.cpp
# g++  -o bin/LN_FAsim -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FAsim.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_NOISEME.o src/LN_NOISEME.cpp
# g++  -o bin/LN_NOISEME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_NOISEME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_GROW_LAYERS.o src/LN_GROW_LAYERS.cpp
# g++  -o bin/LN_GROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GROW_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_3DGROW_LAYERS.o src/LN_3DGROW_LAYERS.cpp
# g++  -o bin/LN_3DGROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_3DGROW_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_DEBUGGING.o src/LN_DEBUGGING.cpp
# g++  -o bin/LN_DEBUGGING -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_DEBUGGING.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_GFACTOR.o src/LN_GFACTOR.cpp
# g++  -o bin/LN_GFACTOR -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GFACTOR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_LEAKY_LAYERS.o src/LN_LEAKY_LAYERS.cpp
# g++  -o bin/LN_LEAKY_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_LEAKY_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_3DCOLUMNS.o src/LN_3DCOLUMNS.cpp
# g++  -o bin/LN_3DCOLUMNS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_3DCOLUMNS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_SHORT_ME.o src/LN_SHORT_ME.cpp
# g++  -o bin/LN_SHORT_ME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SHORT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_FIX_RIM.o src/LN_FIX_RIM.cpp
# g++  -o bin/LN_FIX_RIM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FIX_RIM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_FLOAT_ME.o src/LN_FLOAT_ME.cpp
# g++  -o bin/LN_FLOAT_ME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FLOAT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_IMAGIRO.o src/LN_IMAGIRO.cpp
# g++  -o bin/LN_IMAGIRO -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_IMAGIRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_DIRECT_SMOOTH.o src/LN_DIRECT_SMOOTH.cpp
# g++  -o bin/LN_DIRECT_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_DIRECT_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_RAGRUG.o src/LN_RAGRUG.cpp
# g++  -o bin/LN_RAGRUG -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_RAGRUG.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_EXTREMETR.o src/LN_EXTREMETR.cpp
# g++  -o bin/LN_EXTREMETR -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_EXTREMETR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_TRIAL.o src/LN_TRIAL.cpp
# g++  -o bin/LN_TRIAL -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_TRIAL.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_ZOOM.o src/LN_ZOOM.cpp
# g++  -o bin/LN_ZOOM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_ZOOM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_SMOOTH_RIM.o src/LN_SMOOTH_RIM.cpp
# g++  -o bin/LN_SMOOTH_RIM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SMOOTH_RIM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_COLUMNAR_DIST.o src/LN_COLUMNAR_DIST.cpp
# g++  -o bin/LN_COLUMNAR_DIST -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_COLUMNAR_DIST.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_INTPRO.o src/LN_INTPRO.cpp
# g++  -o bin/LN_INTPRO -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_INTPRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_TEMPSMOOTH.o src/LN_TEMPSMOOTH.cpp
# g++  -o bin/LN_TEMPSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_TEMPSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
# g++    -c -o obj/LN_PHYSIO_PARS.o src/LN_PHYSIO_PARS.cpp
# g++  -o bin/LN_PHYSIO_PARS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_PHYSIO_PARS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
