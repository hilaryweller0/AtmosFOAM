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

volScalarField dbdz
(
    IOobject
    (
        "dbdz",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
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

const dimensionedVector Coriolisf("Coriolisf", dimless/dimTime, viscosity);
//const quaternion CoriolisF(Coriolisf.value());
//const tensor invF(inv(sphericalTensor::I + dt.value()*CoriolisF.R()));

//const scalar Prandtlt = readScalar(laminarTransport.lookup("Prandtlt"));

#include "createFvModels.H"
#include "createFvConstraints.H"
