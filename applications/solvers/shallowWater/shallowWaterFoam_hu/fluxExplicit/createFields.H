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

Info << "Reading in velocity field Uf (velocity on faces)\n" << endl; 
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading in velocity field U (cell centre velocity)\n" << endl; 
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Volume flux through faces
surfaceScalarField volFlux
(
    IOobject("volFlux", runTime.timeName(), mesh),
    Uf & mesh.Sf()
);

surfaceScalarField hup = fvc::interpolate(h, "hup");
surfaceScalarField hc = fvc::interpolate(h, "hc");

// Mass flux through faces
surfaceScalarField flux
(
    IOobject("flux", runTime.timeName(), mesh),
    hup*volFlux
);
flux.oldTime();

