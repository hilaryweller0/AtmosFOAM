Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field phi\n" << endl;

surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

// Read U if present, otherwise create and write (not used)
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

// Read Uf if present, otherwise create and write (not used)
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
    linearInterpolate(fvc::reconstruct(phi))
);
Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
Uf.write();

// Anti-diffusive flux and velocity
surfaceScalarField antiD
(
    IOobject("antiD", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    phi*0
);

surfaceVectorField antiDv
(
    IOobject("antiDv", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    Uf*0
);

