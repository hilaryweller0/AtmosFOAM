Info << "Creating wind field Uf" << endl;
surfaceVectorField Uf
(
    IOobject("Uf", runTime.timeName(), mesh),
    mesh,
    dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
   "fixedValue"
);
