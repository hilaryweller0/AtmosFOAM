#!/bin/bash -e

for dir in linearUpwind* MPDATA* vanLeer* upwind*; do
    ./runOne.sh $dir
done

# Plots
time=600
for dir in linearUpwind* MPDATA* vanLeer* upwind*; do
    ./plotOne.sh $dir $time
done

