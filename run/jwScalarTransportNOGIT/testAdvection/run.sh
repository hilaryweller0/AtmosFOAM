

blockMesh
#setFields
setInitialVelocityField
setAnalyticTracerField -constant
#cp 0/T_analytic 0/T
setVelocityField

#scalarTransportFoamCN 
jwscalarTransportFoam -FE
#paraFoam

#-help can be used to see the options
