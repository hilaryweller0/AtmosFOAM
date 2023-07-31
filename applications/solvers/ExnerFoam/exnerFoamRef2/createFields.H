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
if (!Boussinesq)
{
    Exnerp.oldTime();
}
const label pRefCell = mesh.solution().subDict("solvers").subDict(Exnerp.name())
    .lookupOrDefault<label>("pRefCell", label(0));

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
const volScalarField dThetaadz("dThetaadz", fvc::grad(thetaa) & kdir);

// Update thermodynamics if needed
if (!Boussinesq)
{
    #include "thermoUpdate.H"
}

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
    IOobject("offCentreAdvection", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, ocAlpha)
);
surfaceScalarField aphi = offCentre*phi;

volScalarField rhoRHS
(
    IOobject("rhoRHS", runTime.timeName(), mesh),
    EulerDdt.fvcDdt(rho)
);
rhoRHS.oldTime();

volVectorField Urhs
(
    IOobject("Urhs", runTime.timeName(), mesh),
    (U - U.oldTime())/runTime.deltaT()
);
Urhs.oldTime();
volScalarField thetapRHS
(
    IOobject("thetapRHS", runTime.timeName(), mesh),
    EulerDdt.fvcDdt(rho, thetap)
);

// Fields needed for implicit gravity and acoustic waves
surfaceScalarField gradnTheta(fvc::snGrad(theta));
volScalarField G
(
    "G",
    fvc::average(-offCentre*runTime.deltaT()*rhof*gradnTheta*gSf/(mesh.magSf()*thetaaf))
);
if (!Boussinesq)
{
    G += fvc::average(-offCentre*runTime.deltaT()*rhof*gradnTheta*Cp*fvc::snGrad(Exnerp));
}
G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
volScalarField Psi("Psi", rho/(RbyCv*(Exnera+Exnerp.oldTime())));

