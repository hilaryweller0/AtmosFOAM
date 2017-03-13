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

dimensionedVector Omega(envProperties.lookup("Omega"));
const dimensionedVector Ug
(
    envProperties.lookupOrDefault<dimensionedVector>
    (
        "Ug",
        dimensionedVector("Ug", dimVelocity, vector::zero)
    )
);
