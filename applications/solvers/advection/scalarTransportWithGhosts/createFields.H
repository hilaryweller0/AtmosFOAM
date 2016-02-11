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

Info<< "Reading fields Uf (compulsory), U and phi (optional)\n" << endl;

surfaceVectorField Uf
(
    IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh, IOobject::NO_READ),
    Uf & mesh.Sf()
);

volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::NO_READ),
    fvc::reconstruct(phi)
);

phi.write();
U.write();

// Create fields on the ghost mesh
volScalarField TGhost
(
    IOobject("TGhost", runTime.timeName(), ghostMesh, IOobject::READ_IF_PRESENT),
    ghostMesh.mapToGhost(T)
);

TGhost.write();

surfaceVectorField UfGhost
(
    IOobject("UfGhost", runTime.timeName(), ghostMesh, IOobject::MUST_READ),
    ghostMesh
);

surfaceScalarField phiGhost
(
    IOobject("phiGhost", runTime.timeName(), ghostMesh, IOobject::NO_READ),
    UfGhost & ghostMesh.Sf()
);

// Store div(phi,T) for time-stepping
volScalarField divPhiTGhost = fvc::div(phiGhost, TGhost);
volScalarField divPhiT = ghostMesh.mapFromGhost(divPhiTGhost);
divPhiT.oldTime();


