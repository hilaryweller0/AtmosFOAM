// Read in prognostic variables
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
rv == 0.5*rt0;
volScalarField rl
(
    IOobject
    (
        "rl", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);
rl == 0.1*rt0;

volScalarField T("T", Exner*theta);
volScalarField p("p", pRef*pow(Exner, 1/kappa));
volScalarField es("es", Pcc*pRef*Foam::exp(TccScale*(T - T0)/(T + Tcc)));
volScalarField condenseRate("condenseRate", -rl/dt);

// Latent heat of vapoursation
volScalarField Lv("Lv", Lv0 - (Cpl - Cpv)*(T - T0));

// The flux, including only the gravity term
surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    H.Hdiag()*gd
);

// Moist theta on the face
surfaceScalarField thetaRho
(
    IOobject("thetaf", runTime.timeName(), mesh),
    fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl), "theta")
);

label groundBC = -1;
label topBC = -1;
forAll(mesh.boundaryMesh(), patchi)
{
    if (mesh.boundaryMesh()[patchi].name() == "top")
    {
        topBC = patchi;
    }
    else if (mesh.boundaryMesh()[patchi].name() == "ground")
    {
        groundBC = patchi;
    }
}
if (groundBC == -1 || topBC == -1)
{
    FatalErrorIn("setExnerBalancedH")
        << " no boundary called top or ground. The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);
}


