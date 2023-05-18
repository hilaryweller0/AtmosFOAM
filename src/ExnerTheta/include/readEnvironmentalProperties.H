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
dimensionedScalar magg
(
    max(mag(g), dimensionedScalar("", dimensionSet(0,1,-2,0,0), SMALL))
);
dimensionedVector ghat = g/magg;
dimensionedVector Omega(envProperties.lookup("Omega"));
const surfaceScalarField gSf("gSf", g & mesh.Sf());
const surfaceScalarField ghatSf("ghatSf", ghat & mesh.Sf());
const surfaceScalarField gUnitNormal("gUnitNormal", ghat & (mesh.Sf() / mesh.magSf()));
// Geostrophic Velocity to maintain
const dimensionedVector Ug
(
    envProperties.lookupOrDefault<dimensionedVector>
    (
        "Ug",
        dimensionedVector("Ug", dimVelocity, vector::zero)
    )
);
// Relatation velocity (relax to this velocity in sponge region)
const dimensionedVector Urelax
(
    envProperties.lookupOrDefault<dimensionedVector>
    (
        "Urelax",
        dimensionedVector("Urelax", dimVelocity, vector::zero)
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

// Gravity as a uniformDimensionedVectorField for consistency with
// other turbulence related things
uniformDimensionedVectorField gvf
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    g
);

// i, j and k directions, with k=-ghat
const vector kdir = -ghat.value();
const vector idir = vector(1,0,0);
const vector jdir = -idir ^ kdir;

