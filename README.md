flameletFoam-2.1.x
==================

Description
==================

This code realizes a steady laminar flamelet approach for turbulent non-premixed combustion.
The solver is based on ''rhoReactingFoam'', i.e. it is pressure-based (PISO), compressible and runs with both LES and RAS turbulence.
 
The theory is mainly taken from the work of N. Peters [1, 2] and is based on the view of a turbulent flame as an ensemble of laminar flamelets.
The calculation of these flamelets is a one-dimensional problem and can be done in a pre-processing step.
Integration using a presumed <math>\beta</math>-Probability Density Function (PDF) accounts for the interaction between turbulent fluctuations and chemistry.
The results of the pre-processing procedure are stored in tables which are accessed during the simulation.
Values of interest, such as species mass fraction or enthalpy, are looked-up and connected to the flow using three parameters - the mixture fraction <math>Z</math>, its variance <math>\tilde{Z''^2}</math> and the scalar dissipation rate  <math>\chi</math>.
In doing so, the expensive solution of chemical mechanisms during run-time can be avoided and the run-time thus reduces significantly.

Installation
==================

* Download the package
* Unpack and go to the ''./src'' folder. Set an environment variable to this folder. Alternatively, add this line to your ''.bashrc''. <bash> export LIB_FLAMELET_SRC=$HOME/($YOUR_PATH$)/flameletFoam/src</bash>
* Execute ./Allwmake
