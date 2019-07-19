.. -*- mode: rst -*-

LAYNII version 1.0.0
Tutorials on layering, layer-smoothing, columnar analysis here: https://layerfmri.com/category/code/
This project is licensed under BSD-3-Clause

.. image:: https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png
    :width: 18px
    :target: https://layerfmri.files.wordpress.com/2018/01/sensory_motor_grid.png
    :alt: example image with layers and columns

    
This is set of standalone layer-fMRI C++ programs that do not have any other dependencies. 


Bob Cox and Rick Reynolds wrote few a nii I/O that I recommend. The original version is in the AFNI sources. 
I collected all the necessary files and adapted them for my taste. All the necessary files are::

    nifti1_io.cpp
    nifti2.h
    nifti2_io.h
    nifti_tool.h
    nifticdf.h
    znzlib.h
    nifti1.h
    nifti1_io.h
    nifti2_io.cpp
    nifti_tool.cpp
    nifticdf.cpp
    znzlib.cpp
    
Using linking those allows you to use nii_datatype and load nii files in your own C++ program with the function


    nifti_image * nim=NULL;
    nim = nifti_image_read(filename, 1);

Example
======

My_nii_read.cpp

It reads in a nii file, accesses the data, manipulates the individual voxels writes out the manipulated data as nii


Usage of My_nii_read.cpp
1.) download the all the files with from github E.g. with the command::

    git clone https://github.com/layerfMRI/laynii
    
2.) go into subfolder::

    cd laynii
    
3.) compile it with::

    make all
    
4.) execute it with::

   ./My_nii_read -input input_example.nii -output output.nii -cutoff 3


For more information see: https://layerfmri.com/2017/11/30/using-a-standalone-nii-i-o-in-c/ 

Comment on GSL
======
Parts of LAYNII depend on GSL, thus you should have it installed:


Linux ::

    sudo apt-get install libgsl0-dev


Mac::


    brew install gsl


