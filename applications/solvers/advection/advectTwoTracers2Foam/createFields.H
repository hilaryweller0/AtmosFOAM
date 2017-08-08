Info<< "Reading field T\n" << endl;

volScalarField q1
(
    IOobject
    (
        "q1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Ensure that rho1.oldTime() is stored
q1.oldTime();

volScalarField q1_temp
(
    IOobject
    (
        "q1_temp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);
// Ensure that rho1.oldTime() is stored
q1_temp.oldTime();

volScalarField q2
(
    IOobject
    (
        "q2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
q2.oldTime();

volScalarField q2_temp
(
    IOobject
    (
        "q2_temp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q2
);
// Ensure that rho1.oldTime() is stored
q2_temp.oldTime();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
rho.oldTime();

volScalarField q1_analytic
(
    IOobject
    (
        "q1_analytic",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField q2_analytic
(
    IOobject
    (
        "q2_analytic",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);

volScalarField S1
(
    IOobject
    (
        "S1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);

volScalarField S2
(
    IOobject
    (
        "S2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);

volScalarField S3
(
    IOobject
    (
        "S3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);

volScalarField S4
(
    IOobject
    (
        "S4",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    q1
);

volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

surfaceScalarField flux
(
    IOobject
    (
        "flux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

surfaceScalarField fluxOld
(
    IOobject
    (
        "fluxOld",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


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
    fvc::reconstruct(phi)
);
U.write();

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
Uf.write();

