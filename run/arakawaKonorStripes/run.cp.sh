#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh

setInitialTracerField
mv 0/T 0/theta.anom
mv 0/Tf 0/thetaf.anom

setTheta -CP

setExnerBalancedH -noInterpolate

mv 0/theta 0/theta.bg
mv 0/thetaf 0/thetaf.bg

sumFields -scale0 1 -scale1 1 0 theta 0 theta.bg 0 theta.anom
sumFields -scale0 1 -scale1 1 0 thetaf 0 thetaf.bg 0 thetaf.anom

cp init_0/Uf 0/Uf
set +e
exnerFoamCP

#for t in [0-9]*
#do
#	sumFields -scale0 1 -scale1 -1 $t theta_diff $t theta 0 theta.bg
#	sumFields -scale0 1 -scale1 -1 $t thetaf_diff $t thetaf 0 thetaf.bg
#done
