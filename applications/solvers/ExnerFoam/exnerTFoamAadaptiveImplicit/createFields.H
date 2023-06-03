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
U.oldTime();

volScalarField u
(
    IOobject("u", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);
u = U & idir;
u.oldTime();
volScalarField w
(
    IOobject("w", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);
w = U & kdir;
w.oldTime();

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

Info<< "Reading T\n" << endl;
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
volScalarField theta("theta", T/Exner);
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

surfaceScalarField offCentre
(
    "offCentreAdvection",
    max
    (
        0.5,
        1 - rhof*mesh.magSf()
            /(runTime.deltaT()*mesh.deltaCoeffs()
                *max(mag(phi), dimensionedScalar(phi.dimensions(), SMALL)))
    )
);
surfaceScalarField phiOld("phiOld", (1-offCentre)*phi.oldTime());
surfaceScalarField phiNew("phiNew", offCentre*phi);
volVectorField Urhs
(
    IOobject("Urhs", runTime.timeName(), mesh),
    rho*fvc::weightedReconstruct
    (
        gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf(), 0.25
    ),
    //mesh,
    //dimensionedVector(dimensionSet(1,-2,-2,0,0), vector::zero),
    U.boundaryField().types()
);
Urhs.oldTime();

