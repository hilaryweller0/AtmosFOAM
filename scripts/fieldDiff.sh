#!/bin/bash
if [[ "$#" -lt 2 ]]; then
    echo "Usage: fieldDiff.sh <time> <field>"
    exit 1
fi
sumFields -scale0 1 -scale1 -1 $1 $2_diff 0 $2 $1 $2
