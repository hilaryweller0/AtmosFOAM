// Read in prognostic variables
Info<< "Reading field Exner\n" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Exner.oldTime();

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable";
partNames[1] = "buoyant";

// Set pointer lists for variables in all partitions. Read in later
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);
partitionedVolScalarField rho
(
    "rho", partNames,
    volScalarField
    (
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, scalar(1))
    ),
    sigma
);
partitionedVolScalarField sigmaRho = rho.fraction();

partitionedVolScalarField theta
(
    "theta", partNames, mesh, runTime.timeName(), sigmaRho
);

for(label ip = 0; ip < partNames.size(); ip++)
{
    rho[ip] = pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa);
}
rho.updateSum();
sigmaRho = rho.fraction();
sigmaRho.updateSum();
theta.updateSum();

partitionedSurfaceScalarField sigmaRhof
(
    "sigmaRhof", partNames, fvc::interpolate(sigmaRho[0], "sigmaRho")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRhof.updateSum();

partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, fvc::interpolate(sigma[0], "sigma")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}
sigmaf.updateSum();

partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaRhof
);

partitionedSurfaceScalarField flux
(
    "flux", partNames, Uf.sum() & mesh.Sf()
);
for(label ip = 0; ip < flux.size(); ip++)
{
    flux[ip] = sigmaRhof[ip]*(Uf[ip] & mesh.Sf());
}

partitionedSurfaceScalarField volFlux = flux.divideBy(sigmaRhof, "volFlux");

partitionedVolVectorField u
(
    "u", partNames, fvc::reconstruct(volFlux[0]), sigmaRho
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}

// Store necessary old times and rates of change
volFlux.storeTime();
theta.storeTime();
sigmaRho.storeTime();
flux.sum().oldTime();

Info << "Initialising the drag between partitions" << endl;
surfaceScalarField cloudRadius
(
    "cloudRadius",
    max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1])
);
surfaceScalarField dragCommon = -dragCoeff/cloudRadius
                         * mag(Uf[0] - Uf[1])*(Uf[0] - Uf[1]) & mesh.Sf();
partitionedSurfaceScalarField drag("drag", partNames, dragCommon);

volScalarField Psi("Psi", sigmaRho.sum()/Exner);
Psi.oldTime();

surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    Cp*fvc::interpolate(theta.sum(), "theta")
);

#include "calculateDrag.H"

