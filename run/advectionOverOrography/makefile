#
# Create and run an OpenFOAM case for tracer advection over orography as in
# Schar, Leuenberger, Fuhrer, Luthi and Girard "A New Terrain-Following
# Vertical Coordinate Formulation for Atmospheric Prediction Models", MWR 2002
#

.DEFAULT_GOAL := 10000/Tcontours.pdf

# create the mesh without a mountain
constant/polyMesh/faces : constant/polyMesh/blockMeshDict
	blockMesh > blockMesh.log

# add the mountain
constant/polyMesh/points : constant/polyMesh/faces constant/add2dMountainDict
	add2dMountain > add2dMountain.log

# create a plot of the mesh and view
constant/mesh.pdf : constant/polyMesh/points constant/gmtDicts/mesh
	plotPatchData -constant mesh > plotMesh.log
	gv constant/mesh.pdf &

# initialise the velocity field and the tracer field
0/phi : constant/velocityFieldDict constant/polyMesh/points
	setHorizontalVelocityField
0/T : constant/setScalarOverOrographyDict constant/polyMesh/points
	setScalarOverOrography > setScalar.log

# advect the tracer
10000/T : 0/T 0/phi constant/polyMesh/points system/fvSchemes system/controlDict
	advectionFoam -leapfrog > log

# calculate analytic solution for all time-steps
10000/Texact : constant/setScalarOverOrographyDict
	setScalarOverOrography -tracerFieldFileName Texact >> setScalar.log

# calculate errors
10000/Tdiff : 10000/T 10000/Texact
	for time in [0-9]*; do sumFields $$time Tdiff $$time T $$time Texact -scale1 -1; done > sumFields.log

# plot final contours and errors
10000/Tcontours.pdf : 10000/T 10000/Tdiff constant/gmtDicts/Tcontours
	plotPatchData -time 5000 Tcontours > plotT.log
	plotPatchData -time 10000 Tcontours >> plotT.log
	evince 5000/Tcontours.pdf &
	evince 10000/Tcontours.pdf &

