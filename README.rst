=======
desisim
=======

.. image:: https://img.shields.io/travis/desihub/desisim.svg
    :target: https://travis-ci.org/desihub/desisim
    :alt: Travis Build Status
.. image:: https://coveralls.io/repos/desihub/desisim/badge.svg?service=github
    :target: https://coveralls.io/github/desihub/desisim
    :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/desisim/badge/?version=latest
    :target: http://desisim.readthedocs.org/en/latest/
    :alt: Documentation Status

Introduction
------------

This package contains scripts and packages for simulating DESI spectra.
For full documentation, please visit `desisim on Read the Docs`_

.. _`desisim on Read the Docs`: http://desisim.readthedocs.org/en/latest/

This branch was created to modify the code written by Jaime Ferero-Romano, which is a python wrapper to the fiberassign code.  The fiberassign code, itself, needs a number of input files, for targets, for standard stars, for sky fibers, in addition to those that describe the instrument (locations of positioners and centers of the tiles).  The survey is described by a list of the tiles to be observed.  In desisim, this list comprises a number of epochs.  After each epoch, the list of targets is updated on the basis of observations during the last epoch.  The complete list of targets including sky fibers and standard stars is the mtl (merged target list) file.  The mtl file does not give the redshifts, since it is to represent the real situation, where redshifts aren't know until the observations are made.  The redshifts are provided in a truth file.

To run the code you need to load some modules:

module load specter
module load desispec
module load desitarget
module laod desisim/master

The code is launched with a script, which has as input the designations of the locations of various input and output files.  For example

./my_quicksurvey --output_dir /project/projectdirs/desi/users/rncahn/ --targets_dir /project/projectdirs/desi/users/rncahn/lite  --fiberassign_exec ~/garrigue/fiberassign/bin/fiberassign --epochs_dir ~/garrigue/fiberassign/test  --template_fiberassign ~/garrigue/fiberassign/test/quicksurvey_template_lite.txt --n_epochs 2

In this instance the targets come from a directory with very reduced numbers of galaxies (../rncahn/lite ).  The executable is  ~/garrigue/fiberassign/bin/fiberassign and the code for this is in desihub/fiberassign.  The most confusing aspect is the 'template_fiberassign'.  This is a file required by the executable and begins something like

Targfile /project/projectdirs/desi/mocks/preliminary/mtl/v2/lite/mtl.fits
SStarsfile /project/projectdirs/desi/mocks/preliminary/mtl/v2/lite/stdstars.fits
SkyFfile  /project/projectdirs/desi/mocks/preliminary/mtl/v2/lite/sky.fits
Secretfile /project/projectdirs/desi/mocks/preliminary/mtl/v2/lite/truth.fits

tileFile /project/projectdirs/desi/software/edison/desimodel/0.4/data/footprint/desi-tiles.par
fibFile /project/projectdirs/desi/software/edison/desimodel/0.4/data/focalplane/fiberpos.txt
outDir /project/projectdirs/desi/users/rncahn/tmp/fiberassign 
surveyFile /project/projectdirs/desi/users/rncahn/tmp/survey_list.txt
PrintAscii false
PrintFits true
diagnose false 

The outDir here is actually consistent with output_dir in the call, because the current code adds the two subdirectories!  PrintFits must be true, as it is here.

Here is an example run:

rncahn@edison10:~/desisim/desisim/bin> module load specter
rncahn@edison10:~/desisim/desisim/bin> module load desispec
rncahn@edison10:~/desisim/desisim/bin> module load desitarget
rncahn@edison10:~/desisim/desisim/bin> module load desisim/master
rncahn@edison10:~/desisim/desisim/bin> ./my_quicksurvey --output_dir /project/projectdirs/desi/users/rncahn/ --targets_dir /project/projectdirs/desi/users/rncahn/lite  --fiberassign_exec ~/garrigue/fiberassign/bin/fiberassign --epochs_dir ~/garrigue/fiberassign/test  --template_fiberassign ~/garrigue/fiberassign/test/quicksurvey_template_lite.txt --n_epochs 2
Epoch number 0
mtl epochs [0, 1]
fiber epochs [0]
 tile_ids [  620   621   622 ..., 26664 20903 20920]
10666 tiles to be included in fiberassign
number of truth targets 593965
number of targets 593965
Fri Jul  1 02:15:10 2016 Starting MTL
number of objects in mtl 593965
Fri Jul  1 02:15:13 2016 Finished MTL
Fri Jul  1 02:15:13 2016 Launching fiberassign
Fri Jul  1 02:15:47 2016 Finished fiberassign
Fri Jul  1 02:15:47 2016 104 tiles to gather in zcat
Fri Jul  1 02:15:55 2016 Finished zcat
Epoch number 1
mtl epochs [1]
fiber epochs [1]
 tile_ids [ 9359  9378  9379 ..., 26664 20903 20920]
5333 tiles to be included in fiberassign
Fri Jul  1 02:15:55 2016 Starting MTL
Fri Jul  1 02:16:00 2016 Finished MTL
Fri Jul  1 02:16:00 2016 Launching fiberassign
Fri Jul  1 02:16:20 2016 Finished fiberassign
Fri Jul  1 02:16:20 2016 112 tiles to gather in zcat
Fri Jul  1 02:16:33 2016 Finished zcat


In this instance we had a very small galaxy sample (0<RA,10, -10<DEC <10), which used only 104+112=216 tiles.  The full run was split into two epochs.  

License
-------

desisim is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
