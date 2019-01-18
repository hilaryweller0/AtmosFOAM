Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

// Read in prognostic variables
Info<< "Reading field Uf\n" << endl;
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Cell centre velocity
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
U = fvc::reconstruct(Uf & mesh.Sf());

Info<< "Reading field Exner\n" << endl;
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

Info<< "Reading field theta\n" << endl;
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
thetaf.oldTime();

// Reading in the sponge layer coefficients (if present)
volScalarField muSponge
(
    volScalarField
    (
        IOobject("muSponge", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    )/dt
);
Info << "Maximum sponge value = " << max(muSponge).value()*dt.value()
     << endl;

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

// fixed flux boundary conditions
wordList phiUBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || utype == "fixedValue"
    )
    {
        phiUBCs[patchi] = "fixedValue";
    }
}

Info<< "Calculating face flux fields volFlux and phi\n" << endl;
surfaceScalarField volFlux
(
    IOobject
    (
       "volFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Uf & mesh.Sf()
);
volFlux.oldTime();
surfaceScalarField phi
(
    IOobject
    (
       "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhof*(Uf & mesh.Sf())
);

//  turbulence fields
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Temporary flux with fixed boundary conditions (without hydrostatic terms)
surfaceScalarField volFluxU
(
    IOobject("volFluxU", runTime.timeName(), mesh),
    volFlux,
    phiUBCs
);

// Rates of change for Crank-Nicholson
surfaceScalarField dvolFluxdt
(
    IOobject("dvolFluxdt", runTime.timeName(), mesh),
    gSf - mesh.magSf()*Cp*thetaf*fvc::snGrad(Exner),
    phiUBCs
);

// fix no flow boundaries
forAll(volFluxU.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = U.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
    )
    {
        volFluxU.boundaryFieldRef()[patchi] == 0;
        dvolFluxdt.boundaryFieldRef()[patchi] == 0;
    }
}

