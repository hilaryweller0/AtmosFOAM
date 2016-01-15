#!/bin/bash -v

if [ "$#" -ne 1 ]
then
   echo usage: post.sh case
   exit
fi

case=$1

rm -f $case/Tcontours*.*p* $case/*/Tcontours*.ps $case/*/Tdiff*.p*

# Sum up all of the errors
sumFields -case $case 0 Tdiff 100 Tdiff 200 Tdiff
sumFields -case $case 0 Tdiff 0 Tdiff 300 Tdiff
sumFields -case $case 0 Tdiff 0 Tdiff 400 Tdiff
sumFields -case $case 0 Tdiff 0 Tdiff 500 Tdiff
sumFields -case $case 0 Tdiff 0 Tdiff 600 Tdiff

# plot the sum of the errors
gmtFoam -case $case -time 0 Tdiff
gmtFoam -case $case -time 0 TdiffBig

# Plot the contours every 100s
gmtFoam -case $case -time 0 TcontoursFirst
for time in 100 200 300 400 500; do
    gmtFoam -case $case -time $time TcontoursMid
done
gmtFoam -case $case -time 600 TcontoursLast

# combine all of the output with Tdiff
cat $case/0/Tdiff.ps $case/0/TcontoursFirst.ps $case/100/TcontoursMid.ps \
    $case/200/TcontoursMid.ps $case/300/TcontoursMid.ps \
    $case/400/TcontoursMid.ps $case/500/TcontoursMid.ps \
    $case/600/TcontoursLast.ps  > $case/Tcontours.ps
ps2eps -f $case/Tcontours.ps
epstopdf $case/Tcontours.eps
gv $case/Tcontours.pdf &

# combine all of the output with TdiffBig
cat $case/0/TdiffBig.ps $case/0/TcontoursFirst.ps $case/100/TcontoursMid.ps \
    $case/200/TcontoursMid.ps $case/300/TcontoursMid.ps \
    $case/400/TcontoursMid.ps $case/500/TcontoursMid.ps \
    $case/600/TcontoursLast.ps  > $case/TcontoursBig.ps
ps2eps -f $case/TcontoursBig.ps
epstopdf $case/TcontoursBig.eps
gv $case/TcontoursBig.pdf &

rm $case/Tcontours*.*ps $case/*/Tcontours*.ps $case/*/Tdiff*.ps

#gmtFoam -case $case -time 1000 Tdiff
#gv $case/1000/Tdiff.pdf &

## Plot the solution at certain times
#times="0 120 250 370 500 620 750 870 1000"
#for time in $times; do
#    gmtFoam -case $case -time $time T
#    gv $case/$time/T.pdf &
#done

