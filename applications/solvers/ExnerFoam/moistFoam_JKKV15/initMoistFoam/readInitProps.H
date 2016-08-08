Info<< "Reading initialProperties\n" << endl;

IOdictionary initDict
(
    IOobject
    (
        "initialProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar thetae0(initDict.lookup("thetae"));
const dimensionedScalar rt0(initDict.lookup("totalWaterRatio"));
const dimensionedScalar qt0 = rt0/(1+rt0);
const label maxIters = readLabel(initDict.lookup("maxIters"));
const scalar underRelax = readScalar(initDict.lookup("underRelax"));
const int BCiters  = readLabel(initDict.lookup("BCiters"));
const scalar BCtol = readScalar(initDict.lookup("BCtol", SMALL));
const scalar boundaryRelaxation = readScalar
(
    initDict.lookup("boundaryRelaxation")
);

const scalar thetaPrime = readScalar(initDict.lookup("thetaPrime"));
const point bubbleCentre(initDict.lookup("bubbleCentre"));
const vector bubbleRadii(initDict.lookup("bubbleRadii"));
