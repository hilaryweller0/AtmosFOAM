#!/bin/bash
set -e
rm -rf [0-9]* processor*
blockMesh
terrainFollowingMesh

setTheta

setExnerBalancedH

cp init_0/Uf 0/Uf
createSpongeLayer

exnerFoamH
