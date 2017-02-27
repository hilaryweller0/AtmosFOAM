Info << "\nReading environmentalProperties" << endl;

IOdictionary envProperties
(
    IOobject
    (
        "environmentalProperties",
        runTime.system(),
        rMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar g(envProperties.lookup("g"));

Info << "Reading in the direction of the rotation vector for the Coriolis, "
     << "OmegaHat" << endl;
const vector OmegaHat(envProperties.lookup("OmegaHat"));

Info << "Reading in the magnitude of the Coriolis at y=0, f0" << endl;
const dimensionedScalar f0(envProperties.lookup("f0"));

Info << "Reading in the variation of Coriolis with y, beta" << endl;
const dimensionedScalar beta(envProperties.lookup("beta"));

// Geostrophic Velocity to maintain
const dimensionedVector Ug
(
    envProperties.lookupOrDefault<dimensionedVector>
    (
        "Ug",
        dimensionedVector("Ug", dimVelocity, vector::zero)
    )
);

surfaceVectorField twoOmegaf
(
    IOobject("twoOmegaf", runTime.timeName(), rMesh),
    (f0 + beta * rMesh.Cf().component(1))*OmegaHat
);
twoOmegaf.write();

surfaceVectorField two_dxOmega
(
    IOobject("two_dxOmega", runTime.timeName(), rMesh),
    H.delta() ^ twoOmegaf
);
two_dxOmega.write();

