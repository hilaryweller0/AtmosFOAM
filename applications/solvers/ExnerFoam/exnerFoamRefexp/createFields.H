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

Info<< "Reading perturbation Exner, Exnerp\n" << endl;
volScalarField Exnerp
(
    IOobject
    (
        "Exnerp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading ambient Exner, Exnera\n" << endl;
const volScalarField Exnera
(
    IOobject("Exnera", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Reading theta perturbations, thetap\n" << endl;
volScalarField thetap
(
    IOobject
    (
        "thetap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading ambient theta, thetaa\n" << endl;
const volScalarField thetaa
(
    IOobject("thetaa", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

// thetaa and thetap on the face
const surfaceScalarField thetaaf("thetaaf", linearInterpolate(thetaa));
surfaceScalarField thetapf
(
    IOobject("thetapf", runTime.timeName(), mesh),
    linearInterpolate(thetap)
);
surfaceScalarField thetaf("thetaf", thetaaf + thetapf);
volScalarField theta("theta", thetaa + thetap);
const volScalarField dThetaadz("dThetaadz", fvc::grad(thetaaf) & kdir);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh,
    IOobject::NO_READ, IOobject::AUTO_WRITE),
    pRef*pow(Exnera, 1/RbyCv)/(R*theta)
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
phi.write();

// Update thermodynamics
#include "thermoUpdate.H"

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

volScalarField dRhodt((rho - rho.oldTime())/runTime.deltaT());
dRhodt.oldTime();

volScalarField dRhoThetadt
(
    (rho*thetap - rho.oldTime()*thetap.oldTime())/runTime.deltaT()
);
dRhoThetadt.oldTime();

volVectorField drhoUdt
(
    (rho*U - rho.oldTime()*U.oldTime())/runTime.deltaT()
);
drhoUdt.oldTime();

surfaceScalarField dPhidt((phi - phi.oldTime())/runTime.deltaT());
dPhidt.oldTime();
