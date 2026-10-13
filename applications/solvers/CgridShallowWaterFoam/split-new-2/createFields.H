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
surfaceScalarField hf = fvc::interpolate(h);

Info << "Reading in velocity field U (cell centre velocity)\n" << endl; 
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
);


// Mass flux through faces
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    hf*fvc::flux(U)
);

volVectorField ghGradh
(
    "ghGradh",
    fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf())
);
ghGradh.oldTime();

Info<< "Creating Coriolis Force" << endl;
const volVectorField F("F", ((2.0*Omega) & gHat)*gHat);

mesh.schemes().setFluxRequired(h.name());
