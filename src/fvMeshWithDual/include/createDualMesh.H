Foam::Info
    << "Create dualMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::fvMeshWithDual dualMesh
(
    Foam::IOobject
    (
        "dualMesh",
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    mesh
);
mesh.setDual(dualMesh);
