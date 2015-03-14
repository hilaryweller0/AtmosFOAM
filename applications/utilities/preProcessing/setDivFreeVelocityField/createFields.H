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

    Info<< "Reading field p\n" << endl;
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

    Info<< "Calculating face flux field phi\n" << endl;

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
        Uf & mesh.Sf()
    );

    // Ensure no flow over walls
    for(label iBC = 0; iBC < mesh.boundary().size(); iBC++)
    {
        if (mesh.boundary()[iBC].type() == "wall")
        {
            phi.boundaryField()[iBC] = 0;
        }
    }

