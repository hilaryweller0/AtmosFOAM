volVectorField U
(
    IOobject("U", runTime.timeName(), mesh),
    mesh,
    dimensionedVector("U", dimVelocity, vector(0,0,0)),
    "fixedValue"
);

surfaceVectorField Uf("Uf", linearInterpolate(U));
surfaceScalarField phi("phi", Uf & mesh.Sf());

