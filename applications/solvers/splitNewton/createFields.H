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

volScalarField h_prime
(
    IOobject
    (
        "h_prime",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", h.dimensions(), 0.0),
    h.boundaryField().types()
);

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

volVectorField U_prime
(
    IOobject
    (
        "U_prime",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", U.dimensions(), vector::zero),
    U.boundaryField().types()
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

Info << "Creating E\n" << endl;
volScalarField E
(
    IOobject
    (
        "E",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5*(h)*magSqr(U) + 0.5*magg*sqr(h)
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
