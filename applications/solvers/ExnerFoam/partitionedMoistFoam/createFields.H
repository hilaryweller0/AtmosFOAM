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

Info << "Reading in air and water phases in stable and buoyant partitions"
     << endl;
wordList phaseNames(2);
phaseNames[0] = "air";
phaseNames[1] = "water";
wordList partNames(2);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

partitionedAtmosphere atmosParts(partNames, phaseNames, mesh, thermoDict);

// Select just the stable part
partition& atmos = atmosParts[0];

perfectGasPhase& air = atmos[0].gas();
//fluidSpecie& water = atmos[1];

Info << "Initialising the Exner pressure and the pressure" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Exner", dimless, scalar(1)),
    "partitionedHydrostaticExner"
);
atmosParts.ExnerFromState(Exner);
Exner.oldTime();
Exner.write();

volScalarField p("p", air.pFromExner(Exner));

// Fields needed in the Helmholtz equation
atmosParts.updateCompressibility(Exner);
atmosParts.Psi().oldTime();

surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("gradPcoeff", dimensionSet(1,-1,-2,0,0), scalar(0))
);
atmosParts.setGradPcoeff(gradPcoeff);

surfaceScalarField rhog
(
    IOobject("rhog", runTime.timeName(), mesh),
    fvc::interpolate(atmosParts.rho())*gSf
);

volScalarField rhoDiff
(
    IOobject
    (
        "rhoDiff", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    atmosParts.rho() - atmosParts.Psi()*Exner
);

