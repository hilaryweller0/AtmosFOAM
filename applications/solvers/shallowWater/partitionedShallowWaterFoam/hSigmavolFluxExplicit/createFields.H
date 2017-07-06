// Read in prognostic variables
Info<< "Reading field h\n" << endl;
volScalarField h
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
h.oldTime();

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<volScalarField> hSigma(nParts);
PtrList<volVectorField> U(nParts);
PtrList<surfaceVectorField> Uf(nParts);
PtrList<surfaceScalarField> flux(nParts);
PtrList<surfaceScalarField> volFlux(nParts);
PtrList<surfaceScalarField> hup(nParts);
PtrList<surfaceScalarField> hc(nParts);

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

    hSigma.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"hSigma",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            h*sigma[ip]
        )
    );
    hSigma[ip].oldTime();

    U.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
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
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U[ip])
        )
    );
    Uf[ip].oldTime();

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

    hup.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"hup",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(hSigma[ip],partNames[ip]+"hup")
        )
    );

    hc.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"hc",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(hSigma[ip],partNames[ip]+"hc")
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
            hup[ip]*volFlux[ip]
        )
    );
    flux[ip].oldTime();
}

