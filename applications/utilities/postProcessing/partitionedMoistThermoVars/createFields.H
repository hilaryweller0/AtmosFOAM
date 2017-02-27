Info << "Reading in Exner pressure" << endl;
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

partitionedAtmosphere atmosParts(partNames, phaseNames, mesh, thermoDict, Exner);

// Select just the stable part
//partition& atmos = atmosParts[0];

//perfectGasPhase& air = atmos[0].gas();
//fluidSpecie& water = atmos[1];

