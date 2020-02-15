[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3514298.svg)](https://doi.org/10.5281/zenodo.3514298)

# LAYNII
<img src="https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png" width=350 align="right" />

This is a package of standalone layer-fMRI C++ programs that depends only on a C++ compiler. The main purpose of this package is to provide layer-analysis software that are not (yet) included in the other major MRI analysis software.

Most used programs (so far) are:
- ``LN_3DGROW_LAYERS`` : To generate layer masks based on CSF and WM border lines.
- ``LN_LAYER_SMOOTH`` : For layer-specific spatial smoothing.
- ``LN_BOCO`` : for BOLD correction in VASO.

Tutorials on layering, layer-smoothing, columnar analysis are [here in layerfmri blog](https://layerfmri.com/category/code/).

**Note:** In order to read and write Nifti (.nii) data, I used code that was originally developed from Bob Cox and Rick Reynolds and adapted it for using here.

## Example

``LN_NOISEME.cpp`` reads in a nii file, accesses the data, manipulates the individual voxels, and writes out the manipulated data as nii.

Usage of ``LN_NOISEME.cpp``:

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
4. Execute it with::
```
./LN_NOISEME -input input_example.nii -output Noised.nii -variance 0.4445
```
5. If you want to use LAYNII from anytwhere in you system, you still need to set the paths::

- On Linux add the following to your ~/.bashrc::
```
export PATH="PathToTheLaynioFolder:$PATH"
```
- On Mac add the following to your ~/.bash_profile or ~/.profile::
```
export PATH="PathToTheLaynioFolder:$PATH"
```
- On Windows, you can set paths as follows:
    - On desktop, right-click the very bottom-left corner of the screen to get the Power User Task Menu.
    - From the Power User Task Menu, click System.
    - In the Settings window, scroll down to the related settings section.
    - Click the System info link.
    - In the System window, click the Advanced system settings link in the left navigation pane.
    - In the System Properties window, click on the Advanced tab
    - Then click the Environment Variables button near the bottom of that tab.
    - In the Environment Variables window, highlight the Path variable in the System variables section and click the Edit button.
    - Add or modify the path lines with the paths you want the computer to access. Each different directory is separated with a semicolon, as shown below.
    - There you can add the path to your LAYNII folder

For more information see [this blog post](https://layerfmri.com/2017/11/30/using-a-standalone-nii-i-o-in-c/).

## Comment on cross-platform compatibility
Since January 2020, all remaining dependencies have been removed and LAYNII can be compiled on Linux, Max, and Windows. All you need is a terminal and a C++ compiler.

1. On Linux `g++` is natively included.

2. On Mac, it will be enabled automatically as part of the `comandline developer tools` as soon as you type `g++`` into the terminal. Alternatively, you can also use Xcode.

3. On Windows, a C++ compiler needs to be installed manually. For example with [cygwin](https://cygwin.com/). I followed the instructions in this [video](https://www.youtube.com/watch?v=DAlS4hF_PbY).


## Comment on GSL
Previous versions of LAYNII depend on GSL. I heard your complaints and removed it.


## Comment on makefile and compiler
Some users seemed to have a compiler installed that does not match the actual CPU architecture of the computer. In those cases it can be easier to compile the programs one by one with g++. Copy-paste the following into your terminal instead in Step 3::
```
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
```

# License
This project is licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
