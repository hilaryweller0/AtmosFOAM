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

