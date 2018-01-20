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
partitionedVolScalarFraction sigma
(
    partNames,
    IOobject("sigma", runTime.timeName(), mesh), 
    mesh
);
partitionedVolScalarField rho
(
    partNames,
    IOobject("rho", runTime.timeName(), mesh), 
    mesh,
    sigma
);
partitionedVolScalarFraction sigmaRho = rho.fraction();

partitionedVolScalarField theta
(
    partNames,
    IOobject("theta", runTime.timeName(), mesh), 
    mesh,
    sigmaRho
);

partitionedSurfaceScalarFraction sigmaRhof
(
    partNames,
    fvc::interpolate(sigmaRho.sum(), "sigmaRho")
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
    sigmaRhof[ip].rename(partNames[ip]+".sigmaRhof");
}
sigmaRhof.updateSum();

partitionedSurfaceVectorField Uf
(
    partNames,
    IOobject("Uf", runTime.timeName(), mesh), 
    mesh,
    sigmaRhof
);

partitionedVolVectorField u
(
    partNames,
    IOobject("u", runTime.timeName(), mesh), 
    mesh,
    sigmaRho
);

partitionedSurfaceScalarFraction flux
(
    partNames,
    Uf.sum() & mesh.Sf()
);
for(label ip = 0; ip < flux.size(); ip++)
{
    flux[ip] = sigmaRhof[ip]*(Uf[ip] & mesh.Sf());
}

partitionedSurfaceScalarField volFlux = flux.field(sigmaRhof);

/*partitionedSurfaceScalarField dVolFluxdt

partitionedSurfaceScalarField drag(nParts);




            pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa)

volScalarField Psi("Psi", rhoSum/Exner);
Psi.oldTime();


surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    Cp*fvc::interpolate(rhoTheta)
);

#include "calculateDrag.H"
*/
