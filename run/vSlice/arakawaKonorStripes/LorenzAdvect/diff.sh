#!/bin/bash
reconstructPar

setAnalyticTracerField

for t in [0-9]*
do
	sumFields -scale0 1 -scale1 -1 $t Tf_diff $t Tf $t Tf_analytic
	sumFields -scale0 1 -scale1 -1 $t T_diff $t T $t T_analytic
done

globalSum T_diff
