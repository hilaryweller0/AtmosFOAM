Info<< "Reading U\n" << endl;
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
// theta on the face
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
typedef ThermophysicalTransportModel
<
    compressible::momentumTransportModel, fluidThermo
> tModel;

autoPtr<tModel> thermophysicalTransport(tModel::New(turbulence(), thermo));

surfaceScalarField alphaEff = fvc::interpolate(thermophysicalTransport->kappaEff());
alphaEff.oldTime();

// Reading in the sponge layer coefficients (if present)
const volScalarField muSponge
(
    volScalarField
    (
        IOobject
        (
            "muSponge",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
       )/runTime.deltaT()
);
const tensor spongeT(envProperties.lookup("spongeOrientation"));

const surfaceScalarField muSpongef
(
    "muSpongef",
    (spongeT*fvc::interpolate(muSponge)) && sqr(mesh.Sf()/mesh.magSf())
);

// Used for implicit gravity waves
surfaceScalarField applySIg
(
    IOobject("applySIg", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("applySIg", dimless, scalar(0))

);
surfaceScalarField G
(
    IOobject("G", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("G", dimless, scalar(1))
);

// Fields needed for implicit gravity and acoustic waves
volScalarField Psi("Psi", rho/(RbyCv*Exner));

