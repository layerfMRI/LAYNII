# LayNii

<img src="https://layerfmri.files.wordpress.com/2020/11/laynii_logo_anim.gif"  width="450" align="right" />



This is a package of standalone layer (functional) magnetic resonance imaging (layer-fMRI) C++ programs that depends only on a C++ compiler. The purpose of this package is to provide layer-analysis software that are not (yet) included in the other major MRI analysis software.

Most essential programs (so far) are:
- `LN2_LAYERS`: To generate equi-distant or equi-volume layers from gray matter segmentation. (Alternative to `LN_GROW_LAYERS` in older versions of LayNii).
- `LN_LAYER_SMOOTH`: For layer-specific spatial smoothing.
- `LN_BOCO`: BOLD correction to attain the VASO contrast with SS-SI-VASO.
- `LN2_MULTILATERATE & LN2_PATCH_FLATTEN`: For flattening cortical chunks (see [an application here](https://doi.org/10.1101/2021.11.25.470023 ))


## Citation

If you use LayNii in your research please cite the following article:

- Huber, L., Poser, B. A., Bandettini, P. A., Arora, K., Wagstyl, K., Cho, S., Goense, J., Nothnagel, N., Morgan, A. T., van den Hurk, J., Mueller A. K., Reynolds, R. C., Glen, D. R., Goebel, R. W., Gulban, O. F. (2021). LayNii: A software suite for layer-fMRI. NeuroImage, 118091. https://doi.org/10.1016/j.neuroimage.2021.118091

In addition, please cite the software version of LayNii by using our Zenodo integration:
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3514297.svg)](https://doi.org/10.5281/zenodo.3514297)


## Installation

### Method 1: [Easiest method] Download the zip file
1. Choose the right version for your computer from our releases page: https://github.com/layerfMRI/LAYNII/releases

2. Unzip the downloaded zip file to a desired location (e.g. `/home/user1/LayNii`)

3. Navigate to the unzipped folder in your terminal and execute a LayNii command:
```bash
cd /home/user1/LayNii
./LN2_LAYERS -h
```

### Method 2: Compile yourself
A detailed descriptions of how to set up LayNii is provided here: [https://layerfmri.com/laynii-setup/](https://layerfmri.com/laynii-setup/)
A brief instruction is also given below.

1. Download the latest release and unzip it or clone the repository with the command:
```bash
git clone --depth 1 https://github.com/layerfMRI/laynii
```

2. Change directory to laynii folder:
```bash
cd laynii
```

3. Compile it with:
```bash
make all
```

**Note-1:** See [this comment on cross-platform compatibility](README_APPENDIX.md).

**Note-2:** See [this comment on makefile and compilers](README_APPENDIX.md).

## Docker image

See https://hub.docker.com/repositories/layerfmri .

### Building the docker image

```
git clone --depth 1 https://github.com/layerfMRI/laynii
cd laynii
docker build . -t laynii:latest
```

### Running the docker image

The template below would map some directories on your comptuer
to the `/output` and `/input` inside the container.

`--user "$(id -u):$(id -g)"` would ensure that the output of this docker run
is not owned by the root user (which is the default in docker).

```bash
OUPUT_DIR=FIXME
INPUT_DIR=FIXME

docker run  -it --rm \
            --user "$(id -u):$(id -g)" \
            -v "${OUPUT_DIR}":/output \
            -v "${INPUT_DIR}":/input \
    laynii:latest \
        laynii_command
```

Example to be run from the root folder of the LAYNII repository.

```bash
INPUT_DIR=${PWD}/test_data

mkdir ${PWD}/output
OUTPUT_DIR=${PWD}/output

docker run -it --rm \
           -v "${OUTPUT_DIR}":/output \
           -v "${INPUT_DIR}":/input \
           --user "$(id -u):$(id -g)" \
    laynii:latest \
        LN_BOCO \
        -Nulled "/input/lo_Nulled_intemp.nii" \
        -BOLD "/input/lo_BOLD_intemp.nii" \
        -trialBOCO 40 -shift \
        -output /output/lo_BOLD
```

## Tutorials & Use Cases

Tutorials on layering, layer-smoothing, columnar analysis are [here in layerfmri blog](https://layerfmri.com/category/code/) and [here in ofgulban's youtube channel](https://youtube.com/playlist?list=PLs_umVHtShfadNm8brOweXHUSmqVDTk4q). Various pipeline script in the context of LayNii see the [LayNii_extras](https://github.com/ofgulban/LayNii_extras) Links to instruction of the specific programs are included in the help output of the respective programs and below"

- [LN2_LAYERS algorithm](https://thingsonthings.org/ln2_layers/)
- [LN2_LAYERS example application](https://layerfmri.com/2020/04/24/equivol/)
- [LN2_DEVEIN](https://layerfmri.com/2020/04/02/devein/)
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
- [LN_SKEW](https://layerfmri.com/2020/04/06/qa/)
- [LN_NOISE_KERNEL](https://layerfmri.com/2020/04/06/qa/)
- [LN_LEAKY_LAYERS](https://layerfmri.com/2020/04/24/equivol/)
- [LN_LOITUMA](https://layerfmri.com/2020/04/24/equivol/)
- [LN_GROW_LAYERS usage](https://layerfmri.com/2018/03/11/quick-layering/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2020/04/24/equivol/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2018/07/19/how-to-convert-any-paper-figure-into-a-layer-profile/)
- [LN_GROW_LAYERS example application](https://layerfmri.com/2018/09/26/columns/)
---
## How to contribute?
If you have any issues when using LayNii, or want to request a new feature, we are happy to see them posted on our [issues page](https://github.com/layerfMRI/LayNii/issues). Please employ this as your preferred method (instead of sending individual emails to the authors), since fellow researchers might have similar issues and suggestions.

### Code Contributors
- Renzo Huber
- Omer Faruk Gulban  
- Remi Gau  
- Alessandra Pizzuti  
- Sebastian Dresbach

## License
LayNii is licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

## Acknowledgments
In order to read and write Nifti data, we have adapted code that was originally developed by the Neuroimaging Informatics Technology Initiative. We thank Bob Cox, Daniel Glen and Rick Reynolds. Between 2019-2022, Renzo Huber was supported by NWO Veni Grant in Maastricht University. Since early 2020, development and maintenance of this project is being actively supported by [Brain Innovation](https://www.brainvoyager.com/) as one of the developers ([Omer Faruk Gulban](https://github.com/ofgulban)) works there.

<img src="visuals/supporting.svg" width=1000 />
