[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3514297.svg)](https://doi.org/10.5281/zenodo.3514297)

# LAYNII
<img src="https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png" width=350 align="right" />

This is a package of standalone layer functional magnetic resonance imaging (layer-fMRI) C++ programs that depends only on a C++ compiler. The purpose of this package is to provide layer-analysis software that are not (yet) included in the other major MRI analysis software.

Most used programs (so far) are:
- `LN_GROW_LAYERS` : To generate layer masks based on CSF and WM border lines. Try `LN2_LAYERS` for our newer equi-distant and equi-volume layering algorithm.
- `LN_LAYER_SMOOTH` : For layer-specific spatial smoothing.
- `LN_BOCO` : for BOLD correction in VASO.

Tutorials on layering, layer-smoothing, columnar analysis are [here in layerfmri blog](https://layerfmri.com/category/code/).

## Installation
1. Download the latest release and unzip it or clone the repository with the command:
```
git clone https://github.com/layerfMRI/laynii
```

2. Change directory to laynii folder:
```
cd laynii
```

3. Compile it with:
```
make all
```

## Usage example
For example `LN_NOISEME.cpp` reads in a nii file, accesses the data, manipulates the individual voxels, and writes out the manipulated data as nii. To use `LN_NOISEME.cpp`, `cd` to LAYNII folder and execute the following command in your commandline:
```
./LN_NOISEME -input input_example.nii -output Noised.nii -variance 0.4445
```

### Using LAYNII from a anywhere in your system
If you want to use LAYNII from anywhere in your system, you still need to set the paths:

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

1. On Linux `g++` is included by default.

2. On Mac, it will be enabled automatically as part of the `commandline developer tools` as soon as you type `g++` into the terminal. Alternatively, you can also use Xcode.

3. On Windows, a C++ compiler needs to be installed manually. For example with [cygwin](https://cygwin.com/). I followed the instructions in this [video](https://www.youtube.com/watch?v=DAlS4hF_PbY).

## Comment on makefile and compiler
Some users seemed to have a compiler installed that does not match the actual CPU architecture of the computer. In those cases it can be easier to compile the programs one by one with g++. Copy-paste the following into your terminal instead in Step 3:

```
c++ -std=c++11 -DHAVE_ZLIB  -o LN_BOCO src/LN_BOCO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_3DCOLUMNS src/LN_3DCOLUMNS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_COLUMNAR_DIST src/LN_COLUMNAR_DIST.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_CORREL2FILES src/LN_CORREL2FILES.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_DIRECT_SMOOTH src/LN_DIRECT_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_GRADSMOOTH src/LN_GRADSMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_ZOOM src/LN_ZOOM.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_FLOAT_ME src/LN_FLOAT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_SHORT_ME src/LN_SHORT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_EXTREMETR src/LN_EXTREMETR.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_GFACTOR src/LN_GFACTOR.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_GROW_LAYERS src/LN_GROW_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_IMAGIRO src/LN_IMAGIRO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_INTPRO src/LN_INTPRO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_LEAKY_LAYERS src/LN_LEAKY_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_NOISEME src/LN_NOISEME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_RAGRUG src/LN_RAGRUG.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_SKEW src/LN_SKEW.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_TEMPSMOOTH src/LN_TEMPSMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_TRIAL src/LN_TRIAL.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_PHYSIO_PARS src/LN_PHYSIO_PARS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_INT_ME src/LN_INT_ME.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_LOITUMA src/LN_LOITUMA.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_NOISE_KERNEL src/LN_NOISE_KERNEL.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN2_DEVEIN src/LN2_DEVEIN.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN2_LAYERS src/LN2_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
```
---
# How to contribute?
If you have any issues when using LAYNII, or want to request a new feature, we are happy to see them posted on our [issues page](https://github.com/layerfMRI/LAYNII/issues). Please employ this as your preferred method (instead of  sending individual emails to the authors), since fellow researchers might have similar issues and suggestions.

# License
LAYNII is licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

# Acknowledgments
In order to read and write Nifti (.nii, .nii.gz) data, we have adapted code that was originally developed from Bob Cox and Rick Reynolds.
