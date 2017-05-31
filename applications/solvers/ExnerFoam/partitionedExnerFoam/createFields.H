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
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<volScalarField> theta(nParts);
PtrList<surfaceVectorField> Uf(nParts);

PtrList<volScalarField> rho(nParts);
PtrList<volScalarField> sigmaRho(nParts);
PtrList<volVectorField> u(nParts);
PtrList<surfaceScalarField> flux(nParts);
PtrList<surfaceScalarField> gradPcoeff(nParts);
PtrList<surfaceScalarField> dFluxdt(nParts);

// Read in or set variables in partitions
for (label ip = 0; ip < nParts; ip++)
{
    sigma.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigma",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    sigma[ip].oldTime();

    theta.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"theta",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    theta[ip].oldTime();
    
    Uf.set
    (
        ip,
        new surfaceVectorField
        (
            IOobject
            (
                partNames[ip]+"Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    
    rho.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"rho",
                runTime.timeName(),
                mesh
            ),
            pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa)
        )
    );
    
    sigmaRho.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigmaRho",
                runTime.timeName(),
                mesh
            ),
            sigma[ip]*rho[ip]
        )
    );
    sigmaRho[ip].oldTime();
    
    u.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"u",
                runTime.timeName(),
                mesh
            ),
            H.reconstruct(Uf[ip] & mesh.Sf()),
            Uf[ip].boundaryField().types()
        )
    );
    
    flux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"flux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(sigmaRho[ip])*(Uf[ip] & mesh.Sf())
        )
    );
    flux[ip].oldTime();
    
    gradPcoeff.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"gradPcoeff",
                runTime.timeName(),
                mesh
            ),
            Cp*fvc::interpolate(sigmaRho[ip])*fvc::interpolate(theta[ip])
        )
    );
    
    dFluxdt.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"dFluxdt",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dFluxdt", dimensionSet(1,0,-2,0,0), scalar(0)),
            "fixedValue"
        )
    );
    dFluxdt[ip].oldTime();
}

// Summed variables
volScalarField rhoSum = sigmaRho[0];
for(label ip = 1; ip < nParts; ip++){ rhoSum += sigmaRho[ip]; }
surfaceScalarField rhof = fvc::interpolate(rhoSum);

volScalarField rhoTheta("", sigmaRho[0]*theta[0]);
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigmaRho[ip]*theta[ip];
};

volScalarField Psi("Psi", rhoSum/Exner);
Psi.oldTime();

surfaceScalarField fluxSum("fluxSum", flux[0]);
for(label ip = 1; ip < nParts; ip++){ fluxSum += flux[ip]; }

surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    gradPcoeff[0]
);
for(label ip = 1; ip < nParts; ip++)
{
    gradPcoeffSum += gradPcoeff[ip];
}

