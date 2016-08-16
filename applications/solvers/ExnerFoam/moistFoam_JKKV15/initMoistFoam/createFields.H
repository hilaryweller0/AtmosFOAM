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

Info<< "Reading field thetaRho\n" << endl;
volScalarField thetaRho
(
    IOobject
    (
        "thetaRho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Water vapour and liquid water mass mixing ratios
volScalarField qv
(
    IOobject
    (
        "qv", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ql
(
    IOobject
    (
        "ql", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

// Related mixing ratios
volScalarField qd("qd", 1 - qv - ql);
volScalarField rv("rv", qv/qd);
volScalarField rl("rl", 0*rv);
volScalarField qvs("qvs", qv);

// Thermodynamics variable kappa as a function of mixing ratios
volScalarField cpml("cpml", Cp + Cpv*rv);
volScalarField Rm("Rm", R + Rv*rv);
volScalarField kappam("kappam", Rm/cpml);

// Primitive thermodynamic variables
volScalarField p("p", pFromExner(Exner, kappam, pRef));
volScalarField T("T", TFromThetaRho(thetaRho, Exner, rv, rl, epsilon));
volScalarField rho("rho", moistGasExnerRho(Exner, thetaRho, kappam, R, pRef));

volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));
volScalarField thetaScale("thetaScale", thetaRho/thetaRho);

// Latent heat of vapoursation
volScalarField Lv("Lv", latentHeat(T, Lv0, T0, Cpl, Cpv));

// The flux, including only the gravity term
surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    H.Hdiag()*gd
);


// Moist theta on the face
surfaceScalarField gradPcoeff = gradPCoeff(cpml, thetaRho, rv, epsilon);

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

