IOdictionary rtDict
(
    IOobject
    (
        "totalMoistureDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary velocityFieldDict
(
    IOobject
    (
        "velocityFieldDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

vector d;
vector d0;
scalar r;
scalar rotationAngle;
dimensionedScalar T0("T0", dimensionSet(0,0,0,1,0,0,0), scalar(1.));
dimensionedScalar es0("es0", dimensionSet(1,-1,-2,0,0,0,0), scalar(611.2));
const vector angularVelocity(velocityFieldDict.lookup("solidBodyRotation"));
const vector rotationCentre(velocityFieldDict.lookup("centreOfRotation"));
const vector halfWidth(rtDict.lookup("halfWidth"));
const vector distributionCentre(rtDict.lookup("centre"));
const scalar distributionMax(rtDict.lookupOrDefault<scalar>("maxMagnitude", scalar(0.05)));





