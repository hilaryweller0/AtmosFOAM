# AtmosFOAM
A collection of OpenFOAM computational fluid dynamics applications and libraries for performing atmospheric experiments.  Includes mesh generators, scalar transport and Navier-Stokes solvers, and post-processing and visualisation tools.

## Installation
First install libgdal-dev using apt-get
Install [OpenFOAM dev](https://github.com/OpenFOAM/OpenFOAM-dev).
Compile all AtmosFOAM applications and libraries using `./Allwmake`


## Testing
Unit tests are run automatically using `Allwmake`.  OpenFOAM test cases must be run manually and require [make-common](https://github.com/hertzsprung/make-common) to be installed.
