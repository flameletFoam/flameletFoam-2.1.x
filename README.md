flameletFoam-2.1.x
==================

## Short Description

This code realizes a steady laminar flamelet approach for turbulent non-premixed combustion.
The solver is based on ''rhoReactingFoam'', i.e. it is pressure-based (PISO), compressible and runs with both LES and RAS turbulence.
 
The theory is mainly taken from the work of N. Peters [1, 2] and is based on the view of a turbulent flame as an ensemble of laminar flamelets.
The calculation of these flamelets is a one-dimensional problem and can be done in a pre-processing step.
Integration using a presumed <math>\beta</math>-Probability Density Function (PDF) accounts for the interaction between turbulent fluctuations and chemistry.
The results of the pre-processing procedure are stored in tables which are accessed during the simulation.
Values of interest, such as species mass fraction or enthalpy, are looked-up and connected to the flow using three parameters - the mixture fraction, its variance and the scalar dissipation rate.
In doing so, the expensive solution of chemical mechanisms during run-time can be avoided and the run-time thus reduces significantly.

More information is available on the Extend-bazaar page:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits

## Installation

This version works for OpenFOAM-2.1

* Prepare a directory on your system, e.g.:  

  `mkdir ~/OpenFOAM/flamletFoam/`

* Download the flameletFoam using git:

  `git pull https://github.com/flameletFoam/flameletFoam-2.1.x/ ~/OpenFOAM/flameletFoam/`

* Set an environment variable to the src folder:

  `export LIB_FLAMELET_SRC=~/OpenFOAM/flameletFoam/flameletFoam-2.1.x/src/`

* Execute `./Allwmake`

## Tutorials

There is a RANS and a LES tutorial available:

  `cd ./tutorials/pilotedDiffusionFlame/ras/`
  `./Allrun`

## Notes

More information on using flameletFoam, in particular a description of the workflow including table-generation, is available on:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits

Please feel free to contact me should you find bugs or have suggestions how to make the code better!
