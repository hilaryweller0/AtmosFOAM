//Read dictionaries for settings.
IOdictionary fluidsDict
(
    IOobject
    (
        "setFluidsDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//Settings for numerics.
bool useOperatorSplitting = mesh.solutionDict().lookupOrDefault<bool>
(
    "useOperatorSplitting", false
);
bool semiImplicit = mesh.solutionDict().lookupOrDefault<bool>
(
    "semiImplicit", false
);

//Settings for physics.
bool useCoriolis = mesh.solutionDict().lookupOrDefault<bool>
(
    "useCoriolis", false
);
bool useDrag = mesh.solutionDict().lookupOrDefault<bool>
(
    "useDrag", false
);
bool useDiffusion = mesh.solutionDict().lookupOrDefault<bool>
(
    "useDiffusion", false
);
bool useMassTransfer = mesh.solutionDict().lookupOrDefault<bool>
(
    "useMassTransfer", true
);
bool useBuoyancy = mesh.solutionDict().lookupOrDefault<bool>
(
    "useBuoyancy", true
);
bool useGravity = mesh.solutionDict().lookupOrDefault<bool>
(
    "useGravity", false
);


const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
const int nUCorr = itsDict.lookupOrDefault<int>("nUCorrectors", 1);
const double offCentre = mesh.solutionDict().lookupOrDefault<double>("offCentre", 0.9);
const double offCentreOpSplit = mesh.solutionDict().lookupOrDefault<double>("offCentreOperatorSplitting", 0.5);

dimensionedScalar minimumMassTransfer("minimumMassTransfer",dimensionSet(0,1,0,0,0,0,0),scalar(0));

//Settings for drag.
const int Nmax = itsDict.lookupOrDefault<int>("maximumBubbles", 10);
const double dragCoefficient = itsDict.lookupOrDefault<double>("dragCoefficient", 0.01);
const double sourceMag = itsDict.lookupOrDefault<double>("sourceMagnitude", 0.05);
const vector yNorm(0,1,0);
dimensionedScalar lengthScale("lengthScale",dimensionSet(0,1,0,0,0,0,0),scalar(0.025));
dimensionedScalar bubbleRadiusMin("bubbleRadiusMin",dimensionSet(0,1,0,0,0,0,0),scalar(0));
dimensionedScalar bubbleRadiusMax("bubbleRadiusMax",dimensionSet(0,1,0,0,0,0,0),scalar(1));

dimensionedScalar buoyancyMagnitude("buoyancyMagnitude",dimensionSet(0,1,-1,0,0,0,0),scalar(10));
dimensionedScalar diffusionCoeff("diffusionCoefficient",dimensionSet(0,2,-1,0,0,0,0),scalar(500));
dimensionedScalar timeScale = 10*dt;
