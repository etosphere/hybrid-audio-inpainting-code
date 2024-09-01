
	INTRODUCING SPAIN (SPARSE AUDIO INPAINTER)

This package contains codes for the three variants of SPAIN:

   (1) A-SPAIN
   (2) S-SPAIN OMP
   (3) S-SPAIN H

An illustrative example of performance of all three can be obtained by
running inpainting_main.m in the subfolder SPAIN. This script performs
reconstruction of a signal with single gap and plots the results. All the
parameters are to be set in this script.

Note that the implementation uses toolboxes LTFAT version 2.3.1
(http://ltfat.github.io/) and Sparsify version 0.5 by Thomas Blumensath.
The latter needs not to be installed since the necessary functions are
also part of LTFAT.

Furthermore, three functions from LTFAT / Sparsify need to be modified for
the implementation of S-SPAIN OMP to work:

   (1) franamp.m ........ to compute the norms of atoms of redundant DFT
                          correctly
   (2) greed_omp.m ...... to suppress writing outputs and ensure
                          compatibility with greed_omp_qr.m
   (3) greed_omp_qr.m ... to take complex conjugate atoms into account

All three functions modified such way are provided in this package. If
MATLAB searches for these functions in LTFAT directory instead of the
SPAIN folder, the functions will have to be overwritten in ltfat/frames
(franamp.m), ltfat/thirdparty/sparsify (greed_omp.m) and
ltfat/thirdparty/sparsify/private (greed_omp_qr.m).