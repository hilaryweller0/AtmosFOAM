#!/bin/bash
set -e

CASE=.

blockMesh -case $CASE
add2dMountain -case $CASE
setVelocityField -case $CASE
setScalarOverOrography -case $CASE
