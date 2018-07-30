#!/bin/bash

time=60000

sumFields -scale0 1 -scale1 -1 $time thetafDiff $time thetaf constant thetafRef
gmtFoam -time $time thetafDiff
gv $time/thetafDiff.pdf &

