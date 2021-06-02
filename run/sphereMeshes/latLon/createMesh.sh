#!/bin/sh
#

# Clear out old stuff
rm -rf [0-9]* constant/polyMesh

# Create the mesh
sphPolarLatLonMesh

# Create obj file
writeMeshObj -patchEdges -patchFaces
