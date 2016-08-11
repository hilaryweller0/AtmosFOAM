Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read in air and water phases
wordList phaseNames(2);
phaseNames[0] = "air";
phaseNames[1] = "water";
atmosphere atmos(phaseNames, runTime.timeName(), mesh, thermoDict);

//// Read in water phases
//fluidSpecie water
//(
//    IOobject("waterVapourRho", runTime.timeName(), mesh,
//             IOobject::MUST_READ, IOobject::AUTO_WRITE),
//    IOobject("waterLiquidFrac", runTime.timeName(), mesh,
//             IOobject::MUST_READ, IOobject::AUTO_WRITE),
//    mesh,
//    thermoDict.subDict("waterDict")
//);

