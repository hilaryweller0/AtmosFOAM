#!/bin/sh
#

# Clear out old stuff
rm -rf [0-9]* constant cubeMesh/constant cubeMesh/[0-9]*

# Choose the resolution
NX=6
sed 's/NX/'$NX'/g' cubeMesh/system/blockMeshDictTmp > cubeMesh/system/blockMeshDict

# Create a cube mesh
blockMesh -case cubeMesh

# Take the tangent of the points so that when they are expanded 
# out onto the sphere, it is an equal angle cubed sphere
tanPoints -case cubeMesh

# Extrude out the inner patch into a cubed sphere mesh
extrudeMesh
