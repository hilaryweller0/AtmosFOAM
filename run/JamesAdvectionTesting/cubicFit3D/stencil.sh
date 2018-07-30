#!/bin/bash
set -e
rm -rf [0-9]* processor* *dat *png
blockMesh
setInitialTracerField
setVelocityField
advectionFoam -heun2 | tee advectionFoam.log
./stencil.post.sh
