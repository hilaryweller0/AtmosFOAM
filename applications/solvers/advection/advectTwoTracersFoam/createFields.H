Info<< "Reading field T\n" << endl;

volScalarField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Ensure that T.oldTime() and T.oldTime.oldTime() are stored
T1.oldTime();
T1.oldTime().oldTime();

volScalarField T2
(
    IOobject
    (
        "T2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Ensure that T.oldTime() and T.oldTime.oldTime() are stored
T2.oldTime();
T2.oldTime().oldTime();


volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector::zero),
    "zeroGradient"
);

surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Uf & mesh.Sf()
);

volScalarField k1
(
    IOobject
    (
        "k1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T1)
);

volScalarField k2
(
    IOobject
    (
        "k2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T1)
);

volScalarField k3
(
    IOobject
    (
        "k3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T1)
);

volScalarField k4
(
    IOobject
    (
        "k4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T1)
);
