// Read in prognostic variables
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

// The flux, including only the gravity term
surfaceScalarField un
(
    IOobject("un", runTime.timeName(), mesh),
    H.Hdiag()*gd
);

// Fields needed in the Helmholtz equation
surfaceScalarField gradPcoeff
(
    "gradPcoeff",
    fvc::interpolate
    (
        air.Cp()*theta*atmos.rhoR()
       /(atmos.sumDensity()*air.R()*atmos.volGas())
    )
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
    FatalErrorIn("initMoistFoam_HW")
        << " no boundary called top or ground. The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);
}

// Create Exner to have a fixed value top boundary condition
wordList ExnerBCs(mesh.boundary().size(), "hydrostaticExner");
ExnerBCs[topBC] = "fixedValue";
volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
    atmos.ExnerFromTheta(theta),
    ExnerBCs
);
Exner.boundaryFieldRef()[topBC] == topBCval;

volScalarField T("T", Exner*theta);
volScalarField p("p", air.pFromExner(Exner));
volScalarField rho("rho", atmos.sumDensity());

// Scalaing to create the bubble
volScalarField thetaScale("thetaScale", theta/theta);

