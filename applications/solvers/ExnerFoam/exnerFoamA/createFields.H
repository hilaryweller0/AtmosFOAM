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

thermo.p() = pRef*pow(Exner, 1/kappa);

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
    fvc::interpolate(theta)
);
thermo.T() = theta*Exner;

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

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
    rhof*(linearInterpolate(U) & mesh.Sf())
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
const tensor spongeT(envProperties.lookup("spongeOrientation"));

const surfaceScalarField muSpongef
(
    "muSpongef",
    (spongeT*fvc::interpolate(muSponge)) && sqr(mesh.Sf()/mesh.magSf())
);
Info << "Maximum sponge value = " << (max(muSpongef)*runTime.deltaT()).value()
     << endl;

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

