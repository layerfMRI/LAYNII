## Comment on makefile and compilers
Some users seemed to have a compiler installed that does not match the actual CPU architecture of the computer. In those cases it can be easier to compile the programs with another compiler one by one with g++ (instead of c++).
Some users seemed to have a compiler installed but do not have make installed. Thus, instead of executing 'make all', just copy-paste the following into your terminal in the LayNii folder.

```bash
c++ -std=c++11 -DHAVE_ZLIB -o LN_BOCO src/LN_BOCO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_3DCOLUMNS src/LN_3DCOLUMNS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_COLUMNAR_DIST src/LN_COLUMNAR_DIST.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_CORREL2FILES src/LN_CORREL2FILES.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_DIRECT_SMOOTH src/LN_DIRECT_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_GRADSMOOTH src/LN_GRADSMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_ZOOM src/LN_ZOOM.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_FLOAT_ME src/LN_FLOAT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_SHORT_ME src/LN_SHORT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_EXTREMETR src/LN_EXTREMETR.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_GFACTOR src/LN_GFACTOR.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_GROW_LAYERS src/LN_GROW_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_IMAGIRO src/LN_IMAGIRO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_INTPRO src/LN_INTPRO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_LEAKY_LAYERS src/LN_LEAKY_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_NOISEME src/LN_NOISEME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_RAGRUG src/LN_RAGRUG.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_SKEW src/LN_SKEW.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_TEMPSMOOTH src/LN_TEMPSMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_TRIAL src/LN_TRIAL.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_PHYSIO_PARS src/LN_PHYSIO_PARS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_INT_ME src/LN_INT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_LOITUMA src/LN_LOITUMA.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_NOISE_KERNEL src/LN_NOISE_KERNEL.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_INFO src/LN_INFO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN_CONLAY src/LN_CONLAY.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_DEVEIN src/LN2_DEVEIN.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_RIMIFY src/LN2_RIMIFY.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_LAYERS src/LN2_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_COLUMNS src/LN2_COLUMNS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_CONNECTED_CLUSTERS src/LN2_CONNECTED_CLUSTERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_MULTILATERATE src/LN2_MULTILATERATE.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_PATCH_FLATTEN src/LN2_PATCH_FLATTEN.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_CHOLMO src/LN2_CHOLMO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_PROFILE src/LN2_PROFILE.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB -o LN2_MASK src/LN2_MASK.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz

```
