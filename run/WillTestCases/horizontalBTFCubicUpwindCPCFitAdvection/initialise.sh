#!/bin/bash
set -e

CASE=.

blockMesh -case $CASE
add2dMountain -case $CASE
rm -rf {0,5000}
mkdir -p {0,5000}
setVelocityField -case $CASE
setScalarOverOrography -case $CASE -time 0,5000
mv 5000/T 5000/T_analytic
