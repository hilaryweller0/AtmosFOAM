#!/bin/sh
#

# Clear out old stuff
rm -rf [0-9]* constant [0-9]*

# Choose the resolution
NX=6
sed 's/NX/'$NX'/g' system/blockMeshDictTmp > system/blockMeshDict

# Create a cube mesh
blockMesh

# Take the tangent of the points so that when they are expanded 
# out onto the sphere, it is an equal angle cubed sphere
tanPoints

# Extrude out the inner patch into a cubed sphere mesh
extrudeMesh

# Create obj file
writeMeshObj -patchEdges -patchFaces
