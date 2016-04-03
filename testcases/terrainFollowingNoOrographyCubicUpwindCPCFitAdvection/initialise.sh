#!/bin/bash
set -e

CASE=.

blockMesh -case $CASE
rm -rf {0,10000}
mkdir -p {0,10000}
setVelocityField -case $CASE
setScalarOverOrography -case $CASE -time 0,10000
mv 10000/T 10000/T_analytic
