# LAYNII

<img src="https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png" width=350 align="right" />

This is a package of standalone layer functional magnetic resonance imaging (layer-fMRI) C++ programs that depends only on a C++ compiler. The purpose of this package is to provide layer-analysis software that are not (yet) included in the other major MRI analysis software.

Most used programs (so far) are:
-  `LN2_LAYERS`: To generate equi-distant or equi-volume layers from gray matter segmentation. (Alternative to `LN_GROW_LAYERS` in older versions of LAYNII).
- `LN_LAYER_SMOOTH` : For layer-specific spatial smoothing.
- `LN_BOCO` : for BOLD correction in VASO.

## Citation

If you use LAYNII in your research please cite the following article:

- Huber, L., Poser, B. A., Bandettini, P. A., Arora, K., Wagstyl, K., Cho, S., Goense, J., Nothnagel, N., Morgan, A. T., van den Hurk, J., Reynolds, R. C., Glen, D. R., Goebel, R. W., Gulban, O. F. (2020). LAYNII: A software suite for layer-fMRI. BioRxiv. <https://doi.org/10.1101/2020.06.12.148080>

In addition, please cite the used software version of LAYNII by using our Zenodo integration:
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3514297.svg)](https://doi.org/10.5281/zenodo.3514297)


## Installation
A detailed descriptions of how to set up LAYNII is provided here: [https://layerfmri.com/laynii-setup/](https://layerfmri.com/laynii-setup/)
A brief instruction is also given below.

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

## Tutorials & use cases

Tutorials on layering, layer-smoothing, columnar analysis are [here in layerfmri blog](https://layerfmri.com/category/code/). Various pipeline script in the context of LAYNII see the [LAYNII_extras](https://github.com/ofgulban/LAYNII_extras) Links to instruction of the specific programs are included in the help output of the respective programs and below"

- [LN2_LAYERS algorithm](https://thingsonthings.org/ln2_layers/)
- [LN2_LAYERS example application](https://layerfmri.com/2020/04/24/equivol/)
- [LN_GROW_LAYERS usage](https://layerfmri.com/2018/03/11/quick-layering/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2020/04/24/equivol/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2018/07/19/how-to-convert-any-paper-figure-into-a-layer-profile/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2018/09/26/columns/)
- [LN_BOCO](https://layerfmri.com/2019/03/22/analysispipeline/)
- [LN_COLUMNAR_DIST](https://layerfmri.com/2018/09/26/columns/)
- [LN_IMAGIRO](https://layerfmri.com/2018/09/26/columns/)
- [LN_LAYER_SMOOTH](https://layerfmri.com/2018/11/03/anatomically-informed-spatial-smoothing/)
- [LN_GRADSMOOTH](https://layerfmri.com/2018/11/03/anatomically-informed-spatial-smoothing/)
- [LN_GRADSMOOTH_ITER](https://layerfmri.com/2018/11/03/anatomically-informed-spatial-smoothing/)
- [LN_DIRECT_SMOOTH](https://layerfmri.com/2018/11/03/anatomically-informed-spatial-smoothing/)
- [LN_TEMPSMOOTH](https://layerfmri.com/2018/11/03/anatomically-informed-spatial-smoothing/)
- [LN_INTPRO](https://layerfmri.com/2019/02/05/intensity-projections-in-laynii/)
- [LN_MP2RAGE_DNOISE](https://layerfmri.com/2019/06/22/mp2rage/)
- [LN2_DEVEIN](https://layerfmri.com/2020/04/02/devein/)
- [LN_SKEW](https://layerfmri.com/2020/04/06/qa/)
- [LN_NOISE_KERNEL](https://layerfmri.com/2020/04/06/qa/)
- [LN_LEAKY_LAYERS](https://layerfmri.com/2020/04/24/equivol/)
- [LN_LOITUMA](https://layerfmri.com/2020/04/24/equivol/)
---
## Comment on cross-platform compatibility
Since May 2020, LAYNII is also distributed as pre-compiled binaries for Linux, macOS, and Windows (x32 and x64).
Since January 2020, all remaining dependencies have been removed. This should allow the user to use pre-compiled binaries of LAYNII for the respective operating system. Altenatively LAYNII should also be compilable on Linux, macOS, and Windows. All you need is a terminal and a C++ compiler.

1. On Linux `g++` is included by default.

2. On Mac, it will be enabled automatically as part of the `command line developer tools` as soon as you type `g++` into the terminal. Alternatively, you can also use Xcode.

3. On Windows, a C++ compiler needs to be installed manually. For example with [cygwin](https://cygwin.com/). I followed the instructions in this [video](https://www.youtube.com/watch?v=DAlS4hF_PbY).

## Comment on makefile and compiler
Some users seemed to have a compiler installed that does not match the actual CPU architecture of the computer. In those cases it can be easier to compile the programs with another compiler one by one with g++ (instead of c++).
Some users seemed to have a compiler installed but do not have make installed. Thus, instead of executing 'make all', just copy-paste the following into your terminal in the LAYNII folder.

```
c++ -std=c++11 -DHAVE_ZLIB  -o LN_BOCO src/LN_BOCO.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_3DCOLUMNS src/LN_3DCOLUMNS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_COLUMNAR_DIST src/LN_COLUMNAR_DIST.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_CORREL2FILES src/LN_CORREL2FILES.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_DIRECT_SMOOTH src/LN_DIRECT_SMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_GRADSMOOTH src/LN_GRADSMOOTH.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
c++ -std=c++11 -DHAVE_ZLIB  -o LN_GRADSMOOTH_ITER src/LN_GRADSMOOTH_ITER.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
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
c++ -std=c++11 -DHAVE_ZLIB  -o LN2_LAYERS src/LN2_LAYERS.cpp dep/nifti2_io.cpp dep/znzlib.cpp dep/laynii_lib.cpp -I./dep  -lm -lz
```
---
## How to contribute?
If you have any issues when using LAYNII, or want to request a new feature, we are happy to see them posted on our [issues page](https://github.com/layerfMRI/LAYNII/issues). Please employ this as your preferred method (instead of sending individual emails to the authors), since fellow researchers might have similar issues and suggestions.

## License
LAYNII is licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

## Acknowledgments
In order to read and write Nifti (.nii, .nii.gz) data, we have adapted code that was originally developed by the Neuroimaging Informatics Technology Initiative. We thank Bob Cox, Daniel Glen and Rick Reynolds. Since early 2020, development and maintenance of this project is being actively supported by [Brain Innovation](https://www.brainvoyager.com/) as one of the developers ([Omer Faruk Gulban](https://github.com/ofgulban)) works there.
