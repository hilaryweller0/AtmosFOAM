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
Exnerp.oldTime();

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
thetapf.oldTime();
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

const volScalarField N2
(
    "N2",
    max
    (
        dThetaadz*magg/thetaa,
        dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
    )
);

const surfaceScalarField Ndt
(
    "Ndt", 
    runTime.deltaT()*Foam::sqrt(max
    (
        maxInterp.interpolate(N2), //-gSfhat*fvc::snGrad(thetaa)/thetaaf,
        dimensionedScalar("", dimensionSet(0,0,-2,0,0), SMALL)
    ))
);

/*const surfaceScalarField Ndtf
(
    "Ndtf",
    runTime.deltaT()*Foam::sqrt(max
    (
        -gSfhat*fvc::snGrad(thetaa)/thetaaf,
        dimensionedScalar("", dimensionSet(0,0,-2,0,0), SMALL)
    ))
);*/

volScalarField CoNum
(
    IOobject("CoNum", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(phi/rhof, runTime.deltaT())
);
surfaceScalarField ImEx
(
    "ImEx",
    0.5*(sign(maxInterp.interpolate(CoNum) - CoLimitExp) + 1)
);

/*const surfaceScalarField offCentrev
(
    "offCentrev",
    max(ocAlpha, 1 - 1/(2*Ndt))
    //max(ocAlpha, 0.5*(1 + sqrt(max(sqr(Ndt) - 4,scalar(0)))/(2*Ndt)))
);*/

surfaceScalarField offCentre
(
    IOobject("offCentreAdvection", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, ocAlpha)
);
surfaceScalarField aphi = offCentre*phi;

const surfaceScalarField& offCentrePG = offCentre;
/*(
    IOobject("offCentrePG", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, ocAlpha)
);*/

const surfaceScalarField& offCentreGW = offCentrePG;
/*(
    IOobject("offCentreGW", runTime.timeName(), mesh),
    min(scalar(1), max(scalar(0), 0.5*Ndt))
);*/

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
volScalarField thetapRHS
(
    IOobject("thetapRHS", runTime.timeName(), mesh),
    EulerDdt.fvcDdt(rho, thetap)
);

// Fields needed for implicit gravity and acoustic waves
volScalarField Psi("Psi", rho/(RbyCv*(Exnera+Exnerp.oldTime())));
