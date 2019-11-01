# AtmosFOAM
A collection of OpenFOAM computational fluid dynamics applications and libraries for performing atmospheric simulations.

## Installation

### From source

* Install [OpenFOAM dev](https://github.com/OpenFOAM/OpenFOAM-dev).
* Ensure [AtmosFOAM-tools](https://github.com/AtmosFOAM/AtmosFOAM-tools/) is installed
* Go to directory
cd $WM_PROJECT_USER_DIR
and download AtmosFOAM using:
git clone https://github.com/AtmosFOAM/AtmosFOAM.git
* Export environment variables `~/.bashrc` file:

       export ATMOSFOAM_TOOLS_SRC=/path/to/AtmosFOAM-tools/src
       export GMTU=/path/to/AtmosFOAM-tools/gmtUser
       export ATMOSFOAM_SRC=/path/to/AtmosFOAM/src

* Compile all AtmosFOAM applications and libraries:
cd AtmosFOAM
./Allwmake

