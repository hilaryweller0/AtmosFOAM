# AtmosFOAM
[![DOI](https://zenodo.org/badge/20760151.svg)](https://zenodo.org/badge/latestdoi/20760151)
[![Build Status](https://travis-ci.org/AtmosFOAM/AtmosFOAM.svg?branch=master)](https://travis-ci.org/AtmosFOAM/AtmosFOAM)

A collection of OpenFOAM computational fluid dynamics applications and libraries for performing atmospheric experiments.  Includes mesh generators, scalar transport and Navier-Stokes solvers, and post-processing and visualisation tools.

## Installation

### From source

* First, `apt-get install libgdal-dev`
* Ensure [AtmosFOAM-tools](https://github.com/AtmosFOAM/AtmosFOAM-tools/) is installed
* Install [OpenFOAM dev](https://github.com/OpenFOAM/OpenFOAM-dev).
* Compile all AtmosFOAM applications and libraries using `./Allwmake`

### Ubuntu 17.10 amd64 binaries

    sudo sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -"
    sudo add-apt-repository http://dl.openfoam.org/ubuntu
    sudo add-apt-repository "http://dl.openfoam.org/ubuntu dev"
    sudo add-apt-repository "http://atmosfoam-apt.s3-website-eu-west-1.amazonaws.com dev"
    sudo apt-get update --allow-insecure-repositories
    sudo apt-get install atmosfoam

## Testing
Unit tests are run automatically using `Allwmake`.
