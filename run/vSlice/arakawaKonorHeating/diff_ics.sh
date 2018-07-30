#!/bin/bash
rm -rf 0

setTheta -CP
setExnerBalancedH -noInterpolate
mv 0/theta 0/theta.cp
mv 0/thetaf 0/thetaf.cp
mv 0/Exner 0/Exner.cp

setTheta
setExnerBalancedH
mv 0/theta 0/theta.l
mv 0/thetaf 0/thetaf.l
mv 0/Exner 0/Exner.l

sumFields -scale0 1 -scale1 -1 0 theta_l_cp_diff 0 theta.l 0 theta.cp
sumFields -scale0 1 -scale1 -1 0 thetaf_l_cp_diff 0 thetaf.l 0 thetaf.cp
sumFields -scale0 1 -scale1 -1 0 Exner_l_cp_diff 0 Exner.l 0 Exner.cp
