#!/bin/bash -e


nxs=(50 100 200 400)
dts=(0.01 0.005 0.0025 0.00125)
# Generate the mesh and initialise all cases
for type in orthogonal nonOrthog; do
    i=0
    for nx in ${nxs[*]} ; do
        dt=${dts[$i]}
        echo $type $nx $dt
        ./runAll/init.sh $type $nx $dt
        let i=$i+1
    done
done


# run all test cases
for case in */[1-9]* ; do
    scalarDeformationFoam -case $case >& $case/log &
done

# Cacluate error measures for all cases
for case in */[1-9]*/dt*; do
    ./runAll/calcErrors.sh $case
done

gmtPlot runAll/plotl2.gmt
gmtPlot runAll/plotlinf.gmt

# Plot a load of results
for case in */100x100/dt* */400x400/dt_2.5 ; do
    ./runAll/plotResults.sh $case
done

# plot final errors as a function of dx
mkdir -p runAll/data
echo -e "#dx error1st\n20 0.1\n200 1" > runAll/data/1stOrder.dat
echo -e "#dx error2nd\n20 0.005\n200 0.5" > runAll/data/2ndOrder.dat
echo -e "#dx error3rd\n20 1e-4\n200 0.1" > runAll/data/3rdOrder.dat
DX=10000
for type in orthogonal nonOrthog; do
    for c in 1 10; do
        echo '#dx l2 linf' > $type/errorNorms_c$c.dat
        for nx in 50 100 200 400; do
            dt=`echo $c $nx | awk '{print $1*100/$2}'`
            let dx=$DX/$nx
            case=$type/${nx}x${nx}/dt_$dt
            echo $type $dx $dt $nx $case
            l2=`tail -1 $case/errorNorms.dat | awk '{print $3}'`
            linf=`tail -1 $case/errorNorms.dat | awk '{print $4}'`
            echo $dx $l2 $linf >> $type/errorNorms_c$c.dat
        done
    done
done

gmtPlot runAll/plotl2Convergence.gmt
gmtPlot runAll/plotlinfConvergence.gmt

