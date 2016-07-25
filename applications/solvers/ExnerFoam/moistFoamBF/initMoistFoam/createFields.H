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

volScalarField p("p", pFromExner(Exner, kappa, pRef));
volScalarField T("T", TFromTheta(theta, Exner));

volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));
volScalarField condenseRate("condenseRate", -rl/dt);
volScalarField thetaScale("thetaScale", theta/theta);
volScalarField rvs("rvs", epsilon*(es/(p-es)));

// Latent heat of vapoursation
volScalarField Lv("Lv", latentHeat(T, Lv0, T0, Cpl, Cpv));

// The flux, including only the gravity term
surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    H.Hdiag()*gd
);

// Pressure gradient coefficient to convert grad(Exner) to grad(p)
surfaceScalarField gradPcoeff = gradPCoeff(theta, rv, rl, Cp, epsilon);

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

