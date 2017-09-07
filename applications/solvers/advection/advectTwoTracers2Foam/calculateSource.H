IOdictionary rtDict
(
    IOobject
    (
        "rtDict",
        mesh.time().system(),
        mesh,
        IOobject::READ_IF_PRESENT,
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
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

IOdictionary dict
(
    IOobject
    (
        "advectionDict",
        mesh.time().system(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

vector d;
vector d0;
scalar r;
scalar rotationAngle;
const vector angularVelocity(velocityFieldDict.lookup("solidBodyRotation"));
const vector rotationCentre(velocityFieldDict.lookup("centreOfRotation"));
const vector halfWidth(rtDict.lookup("halfWidth"));
const vector distributionCentre(rtDict.lookup("centre"));
const scalar distributionMax(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(0.05)));

const dimensionedScalar tempScale(dict.lookup("tempScale"));
const dimensionedScalar velocityScale(dict.lookup("velocityScale"));
const dimensionedScalar timeScale(dict.lookup("timeScale"));
const dimensionedScalar es0(dict.lookup("es0"));
const dimensionedScalar muv(dict.lookup("muv"));
const dimensionedScalar mud(dict.lookup("mud"));
const dimensionedScalar Lv(dict.lookup("Lv"));
const dimensionedScalar cp(dict.lookup("cp"));
const dimensionedScalar Rv(dict.lookup("Rv"));
const dimensionedScalar rhoAir(dict.lookup("rhoAir"));

volScalarField es
(
    "es",
    es0*Foam::exp( 17.67*(T-273.15*tempScale)/(T-29.65*tempScale) )
);

const volScalarField rvs
(
    "rvs",
    muv/mud * es/(P-es)
);

