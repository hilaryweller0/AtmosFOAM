Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Creating face flux\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(mesh.Sf().dimensions()*U.dimensions(), 0)
);

volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

dimensionedVector geostrophicWind
(
    "geostrophicWind",
    dimVelocity,
    viscosity
);

dimensionedVector Coriolisf
(
    "Coriolisf",
    dimensionSet(0, 0, -1, 0, 0),
    viscosity
);

const scalar Prandtl = readScalar(viscosity->lookup("Prandtl"));

#include "createFvModels.H"
#include "createFvConstraints.H"
