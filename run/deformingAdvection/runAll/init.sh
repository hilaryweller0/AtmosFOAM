#!/bin/bash -e

if [ "$#" -ne 3 ]
then
   echo usage: init.sh orthogonal'|'nonOrthog nx deltaT
   exit
fi

type=$1
let nx=$2*2
nz=$2
dt=$3
case=$type/${nx}x${nz}

rm -rf $case
mkdir $case

# derived properties
# number of grid points in half of the domain
let NX=$nx/2
let NZ=$NX/2
SQUEEZE=1
EXPAND=1
AM=0
AP=0
if [ "$type" == "orthogonal" ]; then
    echo orthogonal case
elif [ "$type" == "nonOrthog" ]; then
    echo non-orthogonal case
    SQUEEZE=0.5
    EXPAND=2
    AM=-0.1
    AP=0.1
else
    echo in init.sh 'type<orthogonal|nonOrthog>' nx deltaT
    echo type must be one of orthogonal or nonOrthog, not $type
    exit
fi

# Generate the case for a mesh on the sphere with the correct time-step
cp -r dummy/system dummy/constant $case
sed -i -e 's/NX/'$NX'/g' -e 's/NZ/'$NZ'/g' -e 's/AM/'$AM'/g' -e 's/AP/'$AP'/g' \
    -e 's/SQUEEZE/'$SQUEEZE'/g' -e 's/EXPAND/'$EXPAND'/g' \
    $case/constant/polyMesh/blockMeshDict
sed -e 's/DELTAT/'$dt'/g' dummy/system/controlDict > $case/system/controlDict

# Generate the mesh on a plane
blockMesh -case $case

# Create the initial conditions
cp -r dummy/0 $case
setGaussians -case $case
gmtFoam -case $case UT
gv $case/0/UT.pdf &

# wrap the mesh around a cylinder
cylinderMesh -case $case
stitchMesh -overwrite -case $case inlet outlet


