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
theta.oldTime();

// Moist theta on the face
surfaceScalarField thetaRho
(
    IOobject("thetaf", runTime.timeName(), mesh),
    fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl), "theta")
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta*(1+rv/epsilon))*pow(Exner, (1-kappa)/kappa)
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
    (Uf & mesh.Sf())
);

Info << "Calculating V mass circulation field\n" << endl;
surfaceScalarField V
(
    IOobject("V", runTime.timeName(), mesh),
    (Uf & H.delta()),
    VuBCs
);
V.oldTime();

// Rates of change for Crank-Nicholson
surfaceScalarField dVdt
(
    IOobject("dVdt", runTime.timeName(), mesh),
    -(H.delta() & (fvc::interpolate(fvc::div(U,u))+2*(Omega^(Uf-Ug))))
    +(gd - H.magd()*Cp*thetaRho*fvc::snGrad(Exner)),
     VuBCs
);
dVdt.oldTime();

volScalarField divu("divu", fvc::div(Uf & mesh.Sf()));

// Other source terms and thermodynamic variables that need to be stored
volScalarField condenseRate("condenseRate", rv/dt*0);
volScalarField thetaSource("thetaSource", theta/dt*0);
thetaSource.oldTime();
volScalarField rvSource("rvSource", rv/dt*0);
rvSource.oldTime();
volScalarField rlSource("rlSource", rl/dt*0);
rlSource.oldTime();
volScalarField exnerSource("exnerSource", Exner/dt*0);
exnerSource.oldTime();

volScalarField p("p", pRef*pow(Exner, 1/kappa));
volScalarField T("T", theta*Exner);

// fix no flow boundaries
forAll(V.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = U.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
    )
    {
        V.boundaryField()[patchi] == 0;
        dVdt.boundaryField()[patchi] == 0;
    }
}

// thetae for diagnostics
volScalarField Lv("Lv", Lv0 - (Cpl - Cpv)*(T - T0));
volScalarField thetae
(
    IOobject
    (
        "thetae", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    T*pow(p/pRef*epsilon/(rv + epsilon), -R/(Cp+Cpl*(rv+rl)))
   *Foam::exp(Lv*rv/((Cp+Cpl*(rv+rl))*T))
);
thetae.write();
