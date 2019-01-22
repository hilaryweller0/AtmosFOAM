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

partitionedSurfaceScalarField flux
(
    "flux", partNames, linearInterpolate(rho[0])*volFlux[0], sigmaf
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    flux[ip] = linearInterpolate(rho[ip])*volFlux[ip];
}
flux.updateSum();

partitionedSurfaceScalarField massFlux
(
    "massFlux", partNames, linearInterpolate(sigmaRho[0])*volFlux[0], sigmaf
);
for(label ip = 1; ip < massFlux.size(); ip++)
{
    massFlux[ip] = linearInterpolate(sigmaRho[ip])*volFlux[ip];
}
massFlux.updateSum();

partitionedSurfaceScalarField sigmaFlux = volFlux.timesSigma();

Info << "Calculating the cell centre velocity, u, in each partition\n" << endl;
partitionedVolVectorField u
(
    "u", partNames, fvc::reconstruct(volFlux[0]), sigmaRho
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}

Info << "Calculating the cell centre momentum, in each partition\n" << endl;
partitionedVolVectorField momentum
(
    "momentum", partNames, sigmaRho[0]*u[0], sigmaRho
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    momentum[ip] = sigmaRho[ip]*u[ip];
}

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
theta.storeTime();
sigmaRho.storeTime();
flux.sum().oldTime();
flux.storeTime();
rho.storeTime();
sigma.storeTime();
sigmaFlux.storeTime();

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


// Temporary mass and heat transfer fields (for calculations)
const dimensionedScalar transferZero
(
    "transferZero", dimless, scalar(0)
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

Info << "Initialising the mass transfer between partitions\n" << endl;
const dimensionedScalar massTransferZero
(
    "massTransferZero", dimensionSet(1,-3,-1,0,0), scalar(0)
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

Info << "Initialising the theta transfer between partitions\n" << endl;
const dimensionedScalar thetaTransferTermZero
(
    "thetaTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
volScalarField thetaTransferTermTmp
(
    IOobject("thetaTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    thetaTransferTermZero
);
PtrList<partitionedVolScalarField> thetaTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    thetaTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"thetaTransferTerm", partNames, thetaTransferTermTmp
        )
    );
}

Info << "Initialising volFlux transfer between partitions\n" << endl;
const dimensionedScalar volFluxTransferTermZero
(
    "volFluxTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
surfaceScalarField volFluxTransferTermTmp
(
    IOobject("volFluxTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    volFluxTransferTermZero
);
PtrList<partitionedSurfaceScalarField> volFluxTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    volFluxTransferTerm.set
    (
        jp, 
        new partitionedSurfaceScalarField
        (
            partNames[jp]+'.'+"volFluxTransferTerm", partNames, volFluxTransferTermTmp
        )
    );
}


Info << "Initialising mass flux transfer between partitions\n" << endl;
const dimensionedScalar massFluxTransferTermZero
(
    "massFluxTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
surfaceScalarField massFluxTransferTermTmp
(
    IOobject("massFluxTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    massFluxTransferTermZero
);
PtrList<partitionedSurfaceScalarField> massFluxTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massFluxTransferTerm.set
    (
        jp, 
        new partitionedSurfaceScalarField
        (
            partNames[jp]+'.'+"massFluxTransferTerm", partNames, massFluxTransferTermTmp
        )
    );
}


Info << "Initialising velocity transfer between partitions\n" << endl;
const dimensionedScalar velocityTransferTermZero
(
    "velocityTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
volScalarField velocityTransferTermTmp
(
    IOobject("velocityTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    velocityTransferTermZero
);
PtrList<partitionedVolScalarField> velocityTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    velocityTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"velocityTransferTerm", partNames, velocityTransferTermTmp
        )
    );
}


Info << "Initialising momentum transfer between partitions\n" << endl;
const dimensionedScalar momentumTransferTermZero
(
    "momentumTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
volScalarField momentumTransferTermTmp
(
    IOobject("momentumTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    momentumTransferTermZero
);
PtrList<partitionedVolScalarField> momentumTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    momentumTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"momentumTransferTerm", partNames, momentumTransferTermTmp
        )
    );
}

/*Info << "Initialising the velocity transfer between partitions\n" << endl;
const dimensionedVector velocityTransferTermZero
(
    "velocityTransferTermZero", dimensionSet(0,-2,-2,0,0), vector(scalar(0),scalar(0),scalar(0))
);
volVectorField velocityTransferTermTmp
(
    IOobject("velocityTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    velocityTransferTermZero
);
PtrList<partitionedVolVectorField> velocityTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    velocityTransferTerm.set
    (
        jp, 
        new partitionedVolVectorField
        (
            partNames[jp]+'.'+"velocityTransferTerm", partNames, velocityTransferTermTmp
        )
    );
}*/

Info << "Initialising the transfer term between partitions\n" << endl;
volScalarField transferTmp
(
    IOobject("transfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    transferZero
);
PtrList<partitionedVolScalarField> transfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    transfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"transfer", partNames, transferTmp
        )
    );
}

surfaceScalarField transferfTmp
(
    IOobject("transferf", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    transferZero
);
PtrList<partitionedSurfaceScalarField> transferf(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    transferf.set
    (
        jp, 
        new partitionedSurfaceScalarField
        (
            partNames[jp]+'.'+"transferf", partNames, transferfTmp
        )
    );
}

/*
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
);*/

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

