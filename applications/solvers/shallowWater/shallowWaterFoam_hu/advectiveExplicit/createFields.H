// Read in prognostic variables

Info<< "Reading height field, h\n" << endl;
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

Info<< "Reading velocity field, Uf\n" << endl;
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    
Info<< "Reading volume flux field, volFlux\n" << endl;
surfaceScalarField volFlux
(
    IOobject
    (
        "volFlux",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    Uf & mesh.Sf()
);
volFlux.oldTime();

Info<< "Reading velocity field, U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
    //fvc::reconstruct(volFlux),
    //Uf.boundaryField().types()
);
