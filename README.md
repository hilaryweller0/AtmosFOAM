# AtmosFOAM
A collection of OpenFOAM computational fluid dynamics applications and libraries for performing atmospheric simulations.

## Installation

### From source

* Install [OpenFOAM 10](https://github.com/OpenFOAM/OpenFOAM-10).
* Go to directory
cd $WM_PROJECT_USER_DIR
and download AtmosFOAM using:
git clone https://github.com/AtmosFOAM/AtmosFOAM.git
* Export environment variables `~/.bashrc` file:

       export GMTU=/path/to/AtmosFOAM/gmtUser
       export ATMOSFOAM_SRC=/path/to/AtmosFOAM/src

* Compile all AtmosFOAM applications and libraries:
cd AtmosFOAM
./Allwmake

