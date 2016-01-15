Info<< "Reading field Tf\n" << endl;

surfaceScalarField Tf
(
    IOobject
    (
        "Tf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField bf
(
    IOobject
    (
        "bf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Tf * gUnitNormal
);
bf.write();
bf.oldTime();

volVectorField b("b", fvc::reconstruct(bf * mesh.magSf()));
b.oldTime();

// theta is read from disk to pick up the boundary patch types
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(), 
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    b & ghat
);
T.oldTime();
T.write();

Info<< "Reading field Uf\n" << endl;

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

Info << "Reading phi if present\n" << endl;
const surfaceScalarField phi
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
