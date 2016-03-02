# AtmosFOAM
A collection of OpenFOAM computational fluid dynamics applications and libraries for performing atmospheric experiments.  Includes mesh generators, scalar transport and Navier-Stokes solvers, and post-processing and visualisation tools.

## Installation
First, install [OpenFOAM 3.0.0](http://www.openfoam.org/download/) or later.
Compile all AtmosFOAM applications and libraries using `./Allwmake`

## Testing
Unit tests are run automatically using `Allwmake`.  OpenFOAM test cases must be run manually and require [make-common](https://github.com/hertzsprung/make-common) to be installed.
