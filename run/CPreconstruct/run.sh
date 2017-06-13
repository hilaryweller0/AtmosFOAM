#!/bin/bash
set -e
rm -rf [0-9]*
mkdir 0
cp init_0/Tf 0/Tf
testCPreconstruction
