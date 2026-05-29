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
surfaceScalarField hf = linearInterpolate(h);

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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //fvc::reconstruct(Uf & mesh.Sf())
);
U = fvc::reconstruct(Uf & mesh.Sf());
// Mass flux through faces
surfaceScalarField flux
(
    IOobject
    (
        "flux",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    hf*(Uf&mesh.Sf()),
    "fixedValue"
);
flux.boundaryFieldRef() == 0;
flux.write();
flux.oldTime();

// Rate of change of flux
surfaceScalarField dhUdt
(
    IOobject("dhUdt", runTime.name(), mesh),
  - hf*((F^Uf) & mesh.Sf())
  - hf*magg*fvc::snGrad(h+h0)*mesh.magSf()
  - fvc::flux(fvc::div(flux,U))
);
dhUdt.oldTime();

mesh.schemes().setFluxRequired(h.name());

surfaceScalarField U_S_old
(
    IOobject
    (
        "U_S_old",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

surfaceScalarField U_S_new
(
    IOobject
    (
        "U_S_new",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);
