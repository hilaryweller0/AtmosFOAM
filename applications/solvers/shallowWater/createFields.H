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

const volVectorField gHat = mesh.C()/mag(mesh.C());

Info<< "Creating Coriolis Force" << endl;
const volVectorField F("F", ((2.0*Omega) & gHat)*gHat);
volScalarField f("f", (2.0*Omega) & gHat);
f.write();

// Rate of change of momentum, inheriting boundary conditions from U
volVectorField dhUdt
(
    IOobject("dhUdt", runTime.name(), mesh),
    h*U/dt,
    U.boundaryField().types()
);
dhUdt = -fvc::div(phi,U) - h*(F ^ U) - magg*h*fvc::grad(h+h0),
dhUdt -= (gHat & dhUdt)*gHat;
dhUdt.correctBoundaryConditions();
dhUdt.oldTime();

surfaceScalarField ghGradh = magg*hf*fvc::snGrad(h+h0)*mesh.magSf();

mesh.schemes().setFluxRequired(h.name());
