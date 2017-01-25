#!/bin/bash
for t in [0-9]*
do
	sumFields -scale0 1 -scale1 -1 $t theta_diff $t theta 0 theta
	sumFields -scale0 1 -scale1 -1 $t thetaf_diff $t thetaf 0 thetaf
#	sumFields -scale0 1 -scale1 -1 $t Exner_diff $t Exner 0 Exner
done

