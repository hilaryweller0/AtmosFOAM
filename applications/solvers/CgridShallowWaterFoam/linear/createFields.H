// Read in prognostic variables
Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
         runTime.name(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
    ),
    mesh
);
h.oldTime();

Info << "Reading in velocity field Uf (velocity on face)\n" << endl; 
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Create cell centre velocity
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(Uf & mesh.Sf())
);

// Volume flux through faces
surfaceScalarField volFlux
(
    IOobject
    (
        "volFlux",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uf&mesh.Sf(),
    "fixedValue"
);
volFlux.boundaryFieldRef() == 0;
volFlux.write();
volFlux.oldTime();

// Rate of change of flux, inheriting boundary conditions from Uf
surfaceScalarField dUdt
(
    IOobject("dUdt", runTime.name(), mesh),
  - ((F^Uf) & mesh.Sf())
  - magg*fvc::snGrad(h+h0)*mesh.magSf()
);
dUdt.oldTime();

// Rate of change of height
volScalarField dhdt
(
    IOobject("dhdt", runTime.name(), mesh),
    -H*fvc::div(volFlux)
);
dhdt.oldTime();
