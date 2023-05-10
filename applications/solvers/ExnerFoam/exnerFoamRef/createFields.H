Info<< "Reading U and w\n" << endl;
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

/*volScalarField w
(
    IOobject
    (
        "w",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/
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

const volScalarField dThetadz("dThetadz", fvc::grad(thetaa) & ghat);

// Total theta on the face
surfaceScalarField thetaf
(
    IOobject("thetaf", runTime.timeName(), mesh),
    linearInterpolate(thetaa + thetap)
);

#include "thermoUpdate.H"

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    thermo.rho()
    //pRef/(R*(thetaa + thetap))*pow(Exnera + Exnerp, 1/RbyCv)
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
    //rhof*((linearInterpolate(U) & mesh.Sf()) - linearInterpolate(w)*ghatSf)
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
Info << "Maximum sponge value = " << (max(muSponge)*runTime.deltaT()).value()
     << endl;

