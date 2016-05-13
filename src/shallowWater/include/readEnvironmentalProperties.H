Info << "\nReading environmentalProperties" << endl;

IOdictionary environmentalProperties
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar g(environmentalProperties.lookup("g"));
const Switch rotating(environmentalProperties.lookup("rotating"));
const dimensionedScalar beta = rotating ?
                    environmentalProperties.lookup("beta") :
                    dimensionedScalar("beta", -dimTime-dimLength, scalar(0));

