#!/bin/bash

if [ "$#" -ne 2 ]
then
   echo usage: logStats.sh case solvedField
   echo Extracts solution stats \(Courant number and solver iterations\) from log file
   exit
fi

case=$1
f=$2

# Get list of times
times=(`grep 'Time = ' $case/log | awk '{if (NF == 3) print $3}'`)

# Create file of Courant number versus time
c=`grep 'Courant Number mean' $case/log | awk '{print $6}'`

echo "#Time maxC" > $case/c.dat
echo -e ${times[*]}\\n$c | \
    awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS $i: $i) } END{ for (i in a) print a[i] }' >> $case/c.dat
echo Courant numbers written to $case/c.dat

# Create file of iterations versus time
noCorr=`grep nOuterCorrectors $case/system/fvSolution | awk '{print $2}' |\
        awk -F';' '{print $1}'`
niCorr=`grep nCorrectors $case/system/fvSolution | awk '{print $2}' \
        | awk -F';' '{print $1}'`
let nCorr=$noCorr*niCorr
nTimes=${#times[*]}

nIter=(`grep 'Solving for '$f $case/log | awk '{print $15}'`)

echo "#Time n${f}Iter" > $case/n${f}Iter.dat
let it=0
while [ "$it" -lt "$nTimes" ]; do
    let ii=0
    allIters=''
    iterSum=0
    while [ "$ii" -lt "$nCorr" ]; do
        let iterSum=$iterSum+${nIter[$nCorr*$it+$ii]}
        allIters=$(echo $allIters " " ${nIter[$nCorr*$it+$ii]})
        let ii=$ii+1
    done
    echo ${times[$it]} $iterSum $allIters >> $case/n${f}Iter.dat
    let it=$it+1
done
echo Number of iterateions written to $case/n${f}Iter.dat
