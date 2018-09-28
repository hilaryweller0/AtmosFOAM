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

Info << "Reading mountain height, h0\n" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("h0", dimLength, scalar(0))
);
Info << "h0 goes from " << min(h0).value() << " to " << max(h0).value()
     << endl;
const Switch withMountain = max(h0).value() - min(h0).value() > SMALL;

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

