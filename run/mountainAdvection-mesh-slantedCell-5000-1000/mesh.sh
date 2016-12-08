#!/bin/bash -e
rm -rf constant/polyMesh constant/mesh.pdf
blockMesh
slantMesh
checkCellVolumes
setSet -constant -noVTK -batch removeTinyCells
subsetMesh -patch ground -overwrite bigCells
rm -f constant/polyMesh/sets/bigCells constant/polyMesh/sets/tinyCells
checkMesh -constant
if [ -e constant/polyMesh/sets/wrongOrientedFaces ] ; then collapseEdges -constant -overwrite ; fi;
checkMesh -constant
if [ -e constant/polyMesh/sets/zeroAreaFaces ] ; then collapseEdges -constant -overwrite -collapseFaceSet zeroAreaFaces ; fi;

# Remove faces between adjacent slanted cells
topoSet
removeFaces prismFaces -overwrite
combinePatchFaces -overwrite 90

# Plot the mesh
gmtFoam mesh
evince constant/mesh.pdf &


