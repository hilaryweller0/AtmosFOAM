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

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

Info << "Setting rho for each partition\n" << endl;
partitionedVolScalarField rho
(
    "rho", partNames,
    volScalarField
    (
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, scalar(1)),
        sigma[0].boundaryField().types()
    ),
    sigma
);
partitionedVolScalarField sigmaRho = rho.timesSigma();

Info << "Reading in theta in each partition\n" << endl;
partitionedVolScalarField theta
(
    "theta", partNames, mesh, runTime.timeName(), sigmaRho
);

Info << "Calculating density and density in each partition and density weighted theta sum\n" << endl;
for(label ip = 0; ip < partNames.size(); ip++)
{
    rho[ip] = pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa);
}
rho.updateSum();
sigmaRho = rho.timesSigma();
sigmaRho.updateSum();
theta.updateSum();

Info << "Initialising sigmaRhof and sigmaf - values on the faces\n" << endl;
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

Info << "Reading in the velocity in each partition\n" << endl;
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaRhof
);

Info << "Calculating the mass flux and volume flux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux", partNames, Uf[0] & mesh.Sf(), sigmaRhof
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    volFlux[ip] = Uf[ip] & mesh.Sf();
}
volFlux.updateSum();

partitionedSurfaceScalarField flux = volFlux.timesSigma();

Info << "Calculating the cell centre velocity, u, in each partition\n" << endl;
partitionedVolVectorField u
(
    "u", partNames, fvc::reconstruct(volFlux[0]), sigmaRho
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
theta.storeTime();
sigmaRho.storeTime();
flux.sum().oldTime();

Info << "Initialising the compressibility\n" << endl;
volScalarField Psi("Psi", sigmaRho.sum()/Exner);
Psi.oldTime();

surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    Cp*fvc::interpolate(theta.sum(), "theta")
);

Info << "Initialising the drag between partitions" << endl;
surfaceScalarField cloudRadius
(
    "cloudRadius",
    max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1])
);
surfaceScalarField dragCommon = -dragCoeff*sigmaRhof.sum()/cloudRadius
                         * mag(Uf[0] - Uf[1]);

#include "calculateDrag.H"

Info << "Initialising the mass transfer between partitions\n" << endl;
// Temporary mass and heat transfer fields (for calculations)
const dimensionedScalar massTransferZero
(
    "massTransferZero", dimensionSet(1,-3,-1,0,0), scalar(0)
);
const dimensionedScalar dimlessZero
(
    "dimlessZero", dimless, scalar(0)
);
const dimensionedScalar heatTransferZero
(
    "heatTransferZero", dimensionSet(0,-2,0,1,0), scalar(0)
);
const dimensionedScalar densityZero
(
    "densityZero", dimensionSet(1,-3,0,0,0), scalar(0)
);
const dimensionedScalar zeroTemp
(
    "zeroTemp", dimensionSet(0,0,0,1,0), scalar(0)
);
volScalarField massTransferTmp
(
    IOobject("massTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    massTransferZero
);
PtrList<partitionedVolScalarField> massTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"massTransfer", partNames, massTransferTmp
        )
    );
}

// Transfer terms for other variables (assuming nParts == 2)
partitionedVolScalarField transfer
(
    "transfer",
    partNames,
    volScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

partitionedSurfaceScalarField transferf
(
    "transferf",
    partNames,
    surfaceScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

partitionedVolScalarField div
(
    "div",
    partNames,
    fvc::div(flux[0])
);
for(label ip= 1; ip < partNames.size(); ip++)
{
    div[ip] = fvc::div(flux[ip]);
}

partitionedVolScalarField divu
(
    "divu",
    partNames,
    fvc::div(volFlux[0]),
    sigmaRho
);

for(label ip = 1; ip < partNames.size(); ip++)
{
    divu[ip] = fvc::div(volFlux[ip]);
}

