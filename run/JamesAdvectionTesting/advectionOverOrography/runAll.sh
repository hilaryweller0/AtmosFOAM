#!/bin/bash -e
#
# Create and run an OpenFOAM case for tracer advection over orography as in
# Schar, Leuenberger, Fuhrer, Luthi and Girard "A New Terrain-Following
# Vertical Coordinate Formulation for Atmospheric Prediction Models", MWR 2002
#

# create the mesh and add a mountain
blockMesh
add2dMountain

# create a plot of the mesh and view
gmtFoam -constant mesh
evince constant/mesh.pdf &

# initialise the velocity field and the tracer field
setVelocityField
setScalarOverOrography

# advect the tracer
advectionFoam -leapfrog

# plot mid and final contours 
gmtFoam -time 5000 Tcontours
gmtFoam -time 10000 Tcontours
evince 5000/Tcontours.pdf &
evince 10000/Tcontours.pdf &

