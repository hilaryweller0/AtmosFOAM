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

// Store div(phi,T) for time-stepping
volScalarField divPhiT = fvc::div(phi, T);
divPhiT.oldTime();

