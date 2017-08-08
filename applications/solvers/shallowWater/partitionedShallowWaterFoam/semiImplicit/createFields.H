// Read in prognostic variables
Info<< "Reading field h\n" << endl;
volScalarField hSum
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// hH the total h from the Helmholtz eqn
volScalarField hH("hH", hSum);
hH.oldTime();

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<volScalarField> h(nParts);
PtrList<surfaceVectorField> Uf(nParts);
PtrList<volVectorField> u(nParts);
PtrList<surfaceScalarField> volFlux(nParts);
PtrList<surfaceScalarField> hf(nParts);

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
    
    h.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"h",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            hSum*sigma[ip]
        )
    );
    h[ip].oldTime();

    hf.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"hf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            linearInterpolate(h[ip])
        )
    );

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

    volFlux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"volFlux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Uf[ip] & mesh.Sf()
        )
    );
    volFlux[ip].oldTime();

    u.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"u",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::reconstruct(volFlux[ip]),
            Uf[ip].boundaryField().types()
        )
    );
}

// Sum of all of the fluxes in all partitions
surfaceScalarField fluxSum("fluxSum", hf[0]*volFlux[0]);
for(label ip = 1; ip < nParts; ip++)
{
    fluxSum += hf[ip]*volFlux[ip];
}

// Sum hf in all partitions
surfaceScalarField hfSum("hfSum", hf[0]);
for(label ip = 1; ip < nParts; ip++)
{
    hfSum += hf[ip];
}

