// Read in prognostic variables

Info<< "Reading field U, cell centred velocity" << endl;
volVectorField U
(
    IOobject
    (
        "U", runTime.timeName(), mesh, IOobject::MUST_READ
    ),
    mesh
);

surfaceVectorField Uf("Uf", fvc::interpolate(U));

Info<< "Reading h" << endl;
volScalarField h
(
    IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

// h interpolated on to faces
surfaceScalarField hf
(
    IOobject("hf", runTime.timeName(), mesh),
    fvc::interpolate(h)
);

Info<< "Initialising the face flux, phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    hf * (Uf & mesh.Sf())
);

