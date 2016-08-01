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
theta.oldTime();

// Water vapour and liquid water mixing ratios
volScalarField rv
(
    IOobject
    (
        "rv", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField rl
(
    IOobject
    (
        "rl", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);
rl.oldTime();
rv.oldTime();

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

Info<< "Creating field rho from eqn of state\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    moistGasExnerRho(Exner, theta, rv, kappa, R, pRef, epsilon)
);
rho.oldTime();
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

// Fields needed in the Helmholtz equation
volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

surfaceScalarField gradPcoeff = gradPCoeff(theta, rv, rl, Cp, epsilon);

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

// Cell centre velocity
volVectorField u
(
    IOobject("u", runTime.timeName(), mesh),
    H.reconstruct(Uf & mesh.Sf()),
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

Info << "Calculating V mass circulation field\n" << endl;
surfaceScalarField V
(
    IOobject("V", runTime.timeName(), mesh),
    rhof*(Uf & H.delta())
//    phiUBCs
);
V.oldTime();

// Temporary flux with fixed boundary conditions (without hydrostatic terms)
surfaceScalarField Vu
(
    IOobject("Vu", runTime.timeName(), mesh),
    V,
    VuBCs
);

// Velocity flux
//surfaceScalarField Uv("Uv", U/rhof);

// Perpendicular component of the momentum
//surfaceVectorField vperp("vperp", Uf - (Uf & H.delta())*H.delta()/sqr(H.magd()));

// Rates of change for Crank-Nicholson
surfaceScalarField dVdt
(
    IOobject("dVdt", runTime.timeName(), mesh),
    -(H.delta() & (fvc::interpolate(fvc::div(U,u))+2*rhof*(Omega^(Uf-Ug))))
    +rhof*(gd - H.magd()*gradPcoeff*fvc::snGrad(Exner)),
     VuBCs
);
dVdt.oldTime();

volScalarField drhodt("drhodt", -fvc::div(U));
drhodt.oldTime();
volScalarField dRhoThetadt("dRhoThetadt", -fvc::div(U, theta));
dRhoThetadt.oldTime();
volScalarField condenseRate("condenseRate", rv/dt*0);
volScalarField drhovdt("drhovdt", -fvc::div(U, rv));
drhovdt.oldTime();
volScalarField drholdt("drhovdt", -fvc::div(U, rl));
drholdt.oldTime();

// Primitive thermodynamic variables
volScalarField p("p", pFromExner(Exner, kappa, pRef));
volScalarField T("T", TFromTheta(theta, Exner));

// fix no flow boundaries
forAll(Vu.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = U.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
    )
    {
        Vu.boundaryFieldRef()[patchi] = 0;
        dVdt.boundaryFieldRef()[patchi] = 0;
    }
}

// thetae for diagnostics
volScalarField Lv("Lv", latentHeat(T, Lv0, T0, Cpl, Cpv));
volScalarField thetae
(
    IOobject
    (
        "thetae", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    thetaeFromPrimitive(T, p, Lv, rv, rl, pRef, epsilon, Cp, Cpl, R)
);
