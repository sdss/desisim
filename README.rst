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

License
-------

desisim is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
