Info << "\nReading environmentalProperties" << endl;

IOdictionary envProperties
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

dimensionedVector g(envProperties.lookup("g"));
dimensionedVector ghat = g/mag(g);
dimensionedVector Omega(envProperties.lookup("Omega"));
const surfaceScalarField gSf("gSf", g & mesh.Sf());
// Geostrophic Velocity to maintain
const dimensionedVector Ug
(
    envProperties.lookupOrDefault<dimensionedVector>
    (
        "Ug",
        dimensionedVector("Ug", dimVelocity, vector::zero)
    )
);

dimensionedScalar diffusionCoeff
(
    envProperties.lookupOrDefault<dimensionedScalar>
    (
        "diffusionCoefficient",
        dimensionedScalar("K", dimensionSet(0,2,-1,0,0), scalar(0))
    )
);
