#!/bin/sh
#

if [ "$#" -ne 2 ]
then
   echo usage: runningMean.sh window inputFile.dat
   exit
fi

window=$1
inFile=$2
outFile=$3

head -1 $inFile

awk -v n=$window '(FNR>1){for(i=1;i<=NF;++i)
                     {s[i] = s[i] - a[FNR%n,i] + $i; a[FNR%n,i]=$i } }
            (FNR >= n+1){ for(i=1;i<=NF;++i)
                         printf "%s" (i==NF?ORS:OFS), s[i]/n }' \
    $inFile
