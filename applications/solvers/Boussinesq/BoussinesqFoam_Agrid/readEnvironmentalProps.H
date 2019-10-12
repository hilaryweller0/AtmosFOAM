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
const dimensionedScalar nu(envProps.lookup("nu"));
const dimensionedScalar alpha(envProps.lookup("alpha"));

// Uniform Cooling
const dimensionedScalar Q
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "Q", dimensionedScalar(dimensionSet(0,1,-3,0,0), scalar(0))
    )
);

