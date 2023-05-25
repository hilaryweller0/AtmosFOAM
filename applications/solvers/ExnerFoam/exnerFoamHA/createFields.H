Info<< "Reading U, u and w\n" << endl;
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

volScalarField u
(
    IOobject("u", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);
u = U & idir;
volScalarField w
(
    IOobject("w", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);
w = U & kdir;

Info<< "Reading Exner\n" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading theta\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Old part of theta change (and theta on the face)
surfaceScalarField thetaf
(
    IOobject("thetaf", runTime.timeName(), mesh),
    linearInterpolate(theta)
);

#include "thermoUpdate.H"

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    thermo.rho()
);
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

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
    fvc::flux(rho*U)
);
phi.oldTime();

//  turbulence fields
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

// Reading in the sponge layer coefficients (if present)
const volScalarField muSponge
(
    IOobject
    (
        "muSponge",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("muSponge", dimless/dimTime, scalar(0))
);

const surfaceScalarField muSpongef
(
    "muSpongef",
    fvc::interpolate(muSponge)*mag(gUnitNormal)
);
