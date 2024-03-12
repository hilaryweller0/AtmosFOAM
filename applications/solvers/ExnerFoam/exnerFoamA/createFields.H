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
    IOobject("u", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);
u = U & idir;
volScalarField w
(
    IOobject("w", runTime.constant(), mesh, IOobject::MUST_READ),
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
Exner.oldTime();

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

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh,
    IOobject::NO_READ, IOobject::AUTO_WRITE),
    pRef*pow(Exner, 1/RbyCv)/(R*theta)
);
rho.write();
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
phi.write();

#include "thermoUpdate.H"

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

volScalarField CoNum
(
    IOobject("CoNum", runTime.timeName(), mesh),
    CourantNo(phi/rhof, runTime.deltaT())
);
localMax<scalar> maxInterp(mesh);
surfaceScalarField ImEx
(
    "ImEx",
    0.5*(sign(maxInterp.interpolate(CoNum) - CoLimitExp) + 1)
);

surfaceScalarField offCentre
(
    IOobject("offCentreAdvection", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, ocAlpha)
);
surfaceScalarField aphi = offCentre*phi;

surfaceScalarField Ndt
(
    "Ndt", 
    runTime.deltaT()*Foam::sqrt(max
    (
        -gSfhat*fvc::snGrad(theta)/thetaf,
        dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
    ))
);
surfaceScalarField offCentrePG
(
    IOobject("offCentrePG", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, ocAlpha)
);

surfaceScalarField offCentreGW = //offCentrePG;
    max(scalar(0), 1 - 1/max(Ndt, SMALL));

volScalarField rhoRHS
(
    IOobject("rhoRHS", runTime.timeName(), mesh),
    EulerDdt.fvcDdt(rho)
);
rhoRHS.oldTime();

volVectorField Urhs
(
    IOobject("Urhs", runTime.timeName(), mesh),
    EulerDdtv.fvcDdt(rho, U)
);
Urhs.oldTime();
volScalarField thetaRHS
(
    IOobject("thetaRHS", runTime.timeName(), mesh),
    EulerDdt.fvcDdt(rho, theta)
);

// Fields needed for implicit gravity and acoustic waves
volScalarField Psi("Psi", rho/(RbyCv*Exner));

