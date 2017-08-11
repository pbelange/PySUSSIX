# PySUSSIX

A Python wrapper for the frequency analysis tool SUSSIX,
cf. https://cds.cern.ch/record/702438/ .

## Installation on Linux
Run the usual make. If you experience problems with
compiling the Fortran code SUSSIX, try to change the
Makefile's compiler flag to
--fcompiler=gfortran
if you have gfortran installed
(instead of using --fcompiler=gnu95).
