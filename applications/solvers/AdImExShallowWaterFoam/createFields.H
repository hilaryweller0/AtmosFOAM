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

Info<< "Creating Coriolis Force" << endl;
const volVectorField F("F", ((2.0*Omega) & gHat)*gHat);

volVectorField ghGradh = fvc::reconstruct
(
    magg*hf*fvc::snGrad(h+h0)*mesh.magSf()
);
ghGradh -= (ghGradh & gHat)*gHat;

// Rate of change of momentum, inheriting U boundary conditions
volVectorField dhUdt
(
    IOobject("dhUdt", runTime.name(), mesh),
    h*U/dt,
    U.boundaryField().types()
);
dhUdt = -fvc::div(phi,U) - h*(F ^ U) - ghGradh,
dhUdt -= (gHat & dhUdt)*gHat;
dhUdt.oldTime();

mesh.schemes().setFluxRequired(h.name());
