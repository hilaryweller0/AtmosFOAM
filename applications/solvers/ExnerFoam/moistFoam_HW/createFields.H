Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info << "Reading in air and water phases" << endl;
wordList phaseNames(2);
phaseNames[0] = "air";
phaseNames[1] = "water";
atmosphere atmos(phaseNames, mesh, thermoDict);

perfectGasPhase& air = atmos[0].gas();
fluidSpecie& water = atmos[1];

Info << "Reading in theta" << endl;
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

Info << "Calculating, density, Exner, pressure and temperature from atmos and theta" << endl;

volScalarField Exner("Exner", atmos.ExnerFromTheta(theta));
Exner.oldTime();
volScalarField T("T", Exner*theta);
volScalarField p("p", air.pFromExner(Exner));
volScalarField rhoTheta("rhoTheta", air.rho()*theta);
rhoTheta.oldTime();
volScalarField dRhoThetadt("dRhoThetadt", 0*rhoTheta/dt);
dRhoThetadt.oldTime();
volScalarField rhoThetaSource("rhoThetaSource", dRhoThetadt);
surfaceScalarField rhof("rhof", fvc::interpolate(air.rho()));

Info << "Reading in velocity field" << endl;
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

// Fields needed in the Helmholtz equation
volScalarField Psi("Psi", air.rho()/Exner);
Psi.oldTime();
surfaceScalarField gradPcoeff
(
    "gradPcoeff",
    fvc::interpolate(air.Cp()*theta/atmos.volAir())
);

// fixed flux boundary conditions
wordList VuBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "fixedValue";
    }
}

volVectorField rhoU
(
    IOobject("rhoU", runTime.timeName(), mesh),
    air.rho()*H.reconstruct(Uf & mesh.Sf()),
    uBCs
);

Info<< "Calculating face flux field U\n" << endl;
surfaceScalarField U
(
    IOobject
    (
       "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhof*(Uf & mesh.Sf())
//    phiUBCs
);
U.oldTime();

Info << "Calculating V mass circulation field\n" << endl;
surfaceScalarField V
(
    IOobject("V", runTime.timeName(), mesh),
    rhof*(Uf & H.delta())
//    phiUBCs
);
V.oldTime();

// Velocity flux
surfaceScalarField un("un", Uf & mesh.Sf());
volScalarField divu("divu", fvc::div(un));

// Temporary flux with fixed boundary conditions (without hydrostatic terms)
surfaceScalarField Vu
(
    IOobject("Vu", runTime.timeName(), mesh),
    V,
    VuBCs
);

// Rates of change for Crank-Nicholson
surfaceScalarField dVdt
(
    IOobject("dVdt", runTime.timeName(), mesh),
    -(H.delta() & (fvc::interpolate(fvc::div(un, rhoU))+2*rhof*(Omega^(Uf-Ug))))
    +rhof*(gd - H.magd()*gradPcoeff*fvc::snGrad(Exner)),
     VuBCs
);
dVdt.oldTime();

// Reading in the sponge layer coefficients (if present)
surfaceScalarField muSponge
(
    surfaceScalarField
    (
        IOobject("muSponge", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    )/dt
);
Info << "Maximum sponge value = " << max(muSponge).value()*dt.value()
     << endl;

