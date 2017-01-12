// Read in prognostic variables

Info<< "Creating Exner_init\n" << endl;
volScalarField Exner_init
(
    IOobject("Exner_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Setting Exner\n" << endl;
volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
    Exner_init
);

Info<< "Creating theta_init\n" << endl;
volScalarField theta_init
(
    IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Setting theta\n" << endl;
volScalarField theta
(
    IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
    theta_init
);

surfaceScalarField thetaf("thetaf", fvc::interpolate(theta));

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);

surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("U", dimensionSet(0,3,-2,0,0), scalar(0))
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
    FatalErrorIn("setIsothermalBalance")
        << " no boundary called top or ground. The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);
}

