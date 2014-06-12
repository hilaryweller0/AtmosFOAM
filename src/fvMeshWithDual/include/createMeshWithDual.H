Foam::Info
    << "Create mesh with dual for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::fvMeshWithDual mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    true,
    fvMeshWithDual::SPHERICALDIST
);

