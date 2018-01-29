#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh
blockMesh

# Create initial conditions
cp -r init0 0

# run test code
testPartionedField
