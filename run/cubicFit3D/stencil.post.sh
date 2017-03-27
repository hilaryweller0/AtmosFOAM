#!/bin/bash
set -e
head -n64 advectionFoam.log | tail -n36 > stencil.global.dat
head -n2 stencil.global.dat > stencil.global.central.dat
tail -n +3 stencil.global.dat > stencil.global.peripheral.dat

head -n102 advectionFoam.log | tail -n36 > stencil.local.dat
head -n2 stencil.local.dat > stencil.local.central.dat
tail -n +3 stencil.local.dat > stencil.local.peripheral.dat

STENCIL_TYPE=global gnuplot -p stencil.plt
STENCIL_TYPE=local gnuplot -p stencil.plt

