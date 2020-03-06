Info << "\nReading environmentalProperties" << endl;

IOdictionary envProps
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

// Viscosity and thermal diffusivity
//const dimensionedScalar nu(envProps.lookup("nu"));
const dimensionedScalar alpha(envProps.lookup("alpha"));

// Coriolis
const dimensionedVector Omega
(
    envProps.lookupOrDefault<dimensionedVector>
    (
        "Omega",
        dimensionedVector(dimless/dimTime, vector(0,0,0))
    )
);

// RCE
const Switch adjustQ(envProps.lookupOrDefault<Switch>("adjustQ", false));
const dimensionedScalar relaxbTime
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "relaxBuoyancyTimeScale", 
        dimensionedScalar("relaxbTime", dimTime, scalar(0))
    )
);
