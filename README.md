[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3514298.svg)](https://doi.org/10.5281/zenodo.3514298)

# LAYNII
<img src="https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png" width=350 align="right" />

This is a package of standalone layer-fMRI C++ programs that depends only on a C++ compiler. The main purpose of this package is to provide layer-analysis software that are not (yet) included in the other major MRI analysis software.

Most used programs (so far) are:
- ``LN_GROW_LAYERS`` : To generate layer masks based on CSF and WM border lines.
- ``LN_LAYER_SMOOTH`` : For layer-specific spatial smoothing.
- ``LN_BOCO`` : for BOLD correction in VASO.

Tutorials on layering, layer-smoothing, columnar analysis are [here in layerfmri blog](https://layerfmri.com/category/code/).

**Note:** In order to read and write Nifti (.nii) data, I used code that was originally developed from Bob Cox and Rick Reynolds and adapted it for using here.

## Installation
1. Download the all the files with from github E.g. with the command::
```
git clone https://github.com/layerfMRI/laynii
```

2. Change directory to laynii folder:
```
cd laynii
```

3. Compile it with::
```
make all
```

## Usage example
For example `LN_NOISEME.cpp` reads in a nii file, accesses the data, manipulates the individual voxels, and writes out the manipulated data as nii. To use `LN_NOISEME.cpp`, `cd` to LAYNII folder and execute the following command in your commandline:
```
./LN_NOISEME -input input_example.nii -output Noised.nii -variance 0.4445
```

### Using LAYNII from a anywhere in your system
If you want to use LAYNII from anytwhere in your system, you still need to set the paths::

### On Linux
Add the following to your `.bashrc`:
```
export PATH="/path/to/LAYNII:$PATH"
```
### On Mac
Add the following to your `.bash_profile` or `.profile`:
```
export PATH="/path/to/LAYNII:$PATH"
```
### On Windows
You can set paths as follows:
1. On desktop, right-click the very bottom-left corner of the screen to get the Power User Task Menu.
2. From the Power User Task Menu, click System.
3. In the Settings window, scroll down to the related settings section.
4. Click the System info link.
5. In the System window, click the Advanced system settings link in the left navigation pane.
6. In the System Properties window, click on the Advanced tab
7. Then click the Environment Variables button near the bottom of that tab.
8. In the Environment Variables window, highlight the Path variable in the System variables section and click the Edit button.
9. Add or modify the path lines with the paths you want the computer to access. Each different directory is separated with a semicolon, as shown below.
10. There you can add the path to your LAYNII folder

For more information see [this blog post](https://layerfmri.com/2017/11/30/using-a-standalone-nii-i-o-in-c/).

## Comment on cross-platform compatibility
Since January 2020, all remaining dependencies have been removed and LAYNII can be compiled on Linux, Max, and Windows. All you need is a terminal and a C++ compiler.

1. On Linux `g++` is natively included.

2. On Mac, it will be enabled automatically as part of the `comandline developer tools` as soon as you type `g++` into the terminal. Alternatively, you can also use Xcode.

3. On Windows, a C++ compiler needs to be installed manually. For example with [cygwin](https://cygwin.com/). I followed the instructions in this [video](https://www.youtube.com/watch?v=DAlS4hF_PbY).

## Comment on GSL
Previous versions of LAYNII depend on GSL. I heard your complaints and removed it.

## Comment on makefile and compiler
Some users seemed to have a compiler installed that does not match the actual CPU architecture of the computer. In those cases it can be easier to compile the programs one by one with g++. Copy-paste the following into your terminal instead in Step 3::

```
g++ -c -std=c++11 -o obj/nifti2_io.o dep/nifti2_io.cpp
g++ -c -std=c++11 -o obj/nifticdf.o dep/nifticdf.cpp
g++ -c -std=c++11 -o obj/znzlib.o dep/znzlib.cpp
g++ -c -std=c++11 -o obj/laynii_lib.o dep/laynii_lib.cpp
g++ -c -std=c++11 -o  obj/LN_BOCO.o src/LN_BOCO.cpp
g++ -o LN_BOCO -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_BOCO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_MP2RAGE_DNOISE.o src/LN_MP2RAGE_DNOISE.cpp
g++ -o LN_MP2RAGE_DNOISE -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_MP2RAGE_DNOISE.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_LAYER_SMOOTH.o src/LN_LAYER_SMOOTH.cpp
g++ -o LN_LAYER_SMOOTH -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_LAYER_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_3DCOLUMNS.o src/LN_3DCOLUMNS.cpp
g++ -o LN_3DCOLUMNS -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_3DCOLUMNS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_COLUMNAR_DIST.o src/LN_COLUMNAR_DIST.cpp
g++ -o LN_COLUMNAR_DIST -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_COLUMNAR_DIST.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_CORREL2FILES.o src/LN_CORREL2FILES.cpp
g++ -o LN_CORREL2FILES -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_CORREL2FILES.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_DIRECT_SMOOTH.o src/LN_DIRECT_SMOOTH.cpp
g++ -o LN_DIRECT_SMOOTH -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_DIRECT_SMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_GRADSMOOTH.o src/LN_GRADSMOOTH.cpp
g++ -o LN_GRADSMOOTH -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_GRADSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_ZOOM.o src/LN_ZOOM.cpp
g++ -o LN_ZOOM -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_ZOOM.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_FLOAT_ME.o src/LN_FLOAT_ME.cpp
g++ -o LN_FLOAT_ME -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_FLOAT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_SHORT_ME.o src/LN_SHORT_ME.cpp
g++ -o LN_SHORT_ME -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_SHORT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_INT_ME.o src/LN_INT_ME.cpp
g++ -o LN_INT_ME -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_INT_ME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_EXTREMETR.o src/LN_EXTREMETR.cpp
g++ -o LN_EXTREMETR -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_EXTREMETR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_GFACTOR.o src/LN_GFACTOR.cpp
g++ -o LN_GFACTOR -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_GFACTOR.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_GROW_LAYERS.o src/LN_GROW_LAYERS.cpp
g++ -o LN_GROW_LAYERS -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_GROW_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_IMAGIRO.o src/LN_IMAGIRO.cpp
g++ -o LN_IMAGIRO -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_IMAGIRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_INTPRO.o src/LN_INTPRO.cpp
g++ -o LN_INTPRO -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_INTPRO.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_LEAKY_LAYERS.o src/LN_LEAKY_LAYERS.cpp
g++ -o LN_LEAKY_LAYERS -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_LEAKY_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_NOISEME.o src/LN_NOISEME.cpp
g++ -o LN_NOISEME -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_NOISEME.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_RAGRUG.o src/LN_RAGRUG.cpp
g++ -o LN_RAGRUG -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_RAGRUG.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_SKEW.o src/LN_SKEW.cpp
g++ -o LN_SKEW -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_SKEW.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_TEMPSMOOTH.o src/LN_TEMPSMOOTH.cpp
g++ -o LN_TEMPSMOOTH -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_TEMPSMOOTH.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_TRIAL.o src/LN_TRIAL.cpp
g++ -o LN_TRIAL -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_TRIAL.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN_PHYSIO_PARS.o src/LN_PHYSIO_PARS.cpp
g++ -o LN_PHYSIO_PARS -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN_PHYSIO_PARS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o 
g++ -c -std=c++11 -o  obj/LN2_LAYERS.o src/LN2_LAYERS.cpp
g++ -o LN2_LAYERS -std=c++11 -Wall -pedantic -DHAVE_ZLIB -I. obj/LN2_LAYERS.o obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
```

# License
This project is licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
