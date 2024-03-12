// Read in prognostic variables
Info<< "Reading pressure, p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading in buoyancy, b\n" << endl;
volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading in cell centre velocity, U\n" << endl;
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

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf()
);

// Solution optimisation
const dictionary& solDict = mesh.solution().subDict("solution");
const Switch implicitAdvection = solDict.lookup("implicitAdvection");
const int nOuterCorr = solDict.lookupOrDefault<int>("nOuterCorrectors", 2);
const int nCorr = solDict.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr =
        solDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
label pRefCell = 0;
scalar pRefValue = 0.0;

// Gravity direction
const uniformDimensionedVectorField kdir
(
    IOobject("kdir", runTime.constant(), mesh),
    dimensionedVector("kdir", dimless, vector(0,0,1))
);

// Bouyancy on the faces in the z direction
surfaceScalarField bf
(
    "bf",
    fvc::interpolate(b)*mesh.Sf().component(2)*dt
);

// Surface scalar fields for hydrostatic boundary conditions
surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    mesh,
    dt
);


