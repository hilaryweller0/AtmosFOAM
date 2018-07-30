#!/bin/bash
reconstructPar

for t in [0-9]*
do
	sumFields -scale0 1 -scale1 -1 $t theta_diff $t theta 0 theta.bg
	sumFields -scale0 1 -scale1 -1 $t thetaf_diff $t thetaf 0 thetaf.bg
done

