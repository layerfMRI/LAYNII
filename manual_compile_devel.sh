#!/bash

# Dependencies
g++    -c -o obj/nifti2_io.o nifti2_io.cpp
g++    -c -o obj/nifticdf.o nifticdf.cpp
g++    -c -o obj/znzlib.o znzlib.cpp
g++    -c -o obj/laynii_lib.o laynii_lib.cpp

# High priority LAYNII programs
g++    -c -o obj/LN_MP2RAGE_DNOISE.o LN_MP2RAGE_DNOISE.cpp -O2
g++  -o LN_MP2RAGE_DNOISE -Wall -pedantic -DHAVE_ZLIB -I. -O2  obj/LN_MP2RAGE_DNOISE.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
g++    -c -o obj/LN_BOCO.o LN_BOCO.cpp
g++  -o LN_BOCO -Wall -pedantic -DHAVE_ZLIB -I. -O2  obj/LN_BOCO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
g++    -c -o obj/LN_LAYER_SMOOTH.o LN_LAYER_SMOOTH.cpp
g++  -o LN_LAYER_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_LAYER_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
g++    -c -o obj/LN_3DGROW_LAYERS.o LN_3DGROW_LAYERS.cpp
g++  -o LN_3DGROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_3DGROW_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o

# Other LAYNII programs
# g++    -c -o obj/LN_CORREL2FILES.o LN_CORREL2FILES.cpp
# g++  -o LN_CORREL2FILES -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_CORREL2FILES.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_DEBUGGING.o LN_DEBUGGING.cpp
# g++  -o LN_DEBUGGING -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_DEBUGGING.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_DIRECT_SMOOTH.o LN_DIRECT_SMOOTH.cpp
# g++  -o LN_DIRECT_SMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_DIRECT_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_EXTREMETR.o LN_EXTREMETR.cpp
# g++  -o LN_EXTREMETR -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_EXTREMETR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_GFACTOR.o LN_GFACTOR.cpp
# g++  -o LN_GFACTOR -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GFACTOR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_GRADSMOOTH.o LN_GRADSMOOTH.cpp
# g++  -o LN_GRADSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GRADSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o

# g++    -c -o obj/LN_GROW_LAYERS.o LN_GROW_LAYERS.cpp
# g++  -o LN_GROW_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_GROW_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_IMAGIRO.o LN_IMAGIRO.cpp
# g++  -o LN_IMAGIRO -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_IMAGIRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_INTPRO.o LN_INTPRO.cpp
# g++  -o LN_INTPRO -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_INTPRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_LEAKY_LAYERS.o LN_LEAKY_LAYERS.cpp
# g++  -o LN_LEAKY_LAYERS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_LEAKY_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_NOISEME.o LN_NOISEME.cpp
# g++  -o LN_NOISEME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_NOISEME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_RAGRUG.o LN_RAGRUG.cpp
# g++  -o LN_RAGRUG -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_RAGRUG.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_SHORT_ME.o LN_SHORT_ME.cpp
# g++  -o LN_SHORT_ME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SHORT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_SKEW.o LN_SKEW.cpp
# g++  -o LN_SKEW -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SKEW.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_FAsim.o LN_FAsim.cpp
# g++    -c -o obj/LN_SMOOTH_RIM.o LN_SMOOTH_RIM.cpp
# g++  -o LN_SMOOTH_RIM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_SMOOTH_RIM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_TEMPSMOOTH.o LN_TEMPSMOOTH.cpp
# g++  -o LN_TEMPSMOOTH -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_TEMPSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_TRIAL.o LN_TRIAL.cpp
# g++  -o LN_TRIAL -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_TRIAL.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o

# g++  -o LN_FAsim -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FAsim.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_3DCOLUMNS.o LN_3DCOLUMNS.cpp
# g++  -o LN_3DCOLUMNS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_3DCOLUMNS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_FIX_RIM.o LN_FIX_RIM.cpp
# g++  -o LN_FIX_RIM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FIX_RIM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_FLOAT_ME.o LN_FLOAT_ME.cpp
# g++  -o LN_FLOAT_ME -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_FLOAT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_COLUMNAR_DIST.o LN_COLUMNAR_DIST.cpp
# g++  -o LN_COLUMNAR_DIST -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_COLUMNAR_DIST.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
# g++    -c -o obj/LN_ZOOM.o LN_ZOOM.cpp
# g++  -o LN_ZOOM -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_ZOOM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o

# g++    -c -o obj/LN_PHYSIO_PARS.o LN_PHYSIO_PARS.cpp
# g++  -o LN_PHYSIO_PARS -Wall -pedantic -DHAVE_ZLIB -I.  obj/LN_PHYSIO_PARS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
