#!/bin/bash
set -e
rm -rf [0-9]* processor*
blockMesh

setInitialTracerField

cp init_0/Uf 0/Uf

#decomposePar -force -constant
#mpirun --hostfile machines -np 2 advectiveFoamF -parallel
advectiveFoamF
./diff.sh
