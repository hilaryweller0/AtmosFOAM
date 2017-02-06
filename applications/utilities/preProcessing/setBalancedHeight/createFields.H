// Read in prognostic variables

Info<< "Reading field Uf, velocity on the faces" << endl;
surfaceVectorField Uf
(
    IOobject
    (
        "Uf", runTime.timeName(), mesh, IOobject::MUST_READ
    ),
    mesh
);

Info<< "Reading h" << endl;
volScalarField h
(
    IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

// h interpolated on to faces
surfaceScalarField hf
(
    IOobject("hf", runTime.timeName(), mesh),
    fvc::interpolate(h)
);

// fixed flux boundary conditions
wordList VuBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    const word htype = h.boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "fixedValue";
    }
}

Info << "Calculating the cell centre velocity, u" << endl;
volVectorField u
(
    IOobject("u", runTime.timeName(), mesh,
             IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
    H.reconstruct(Uf & mesh.Sf()),
    uBCs
);
u.write();

Info<< "Initialising the absolute face flux, U" << endl;
surfaceScalarField U
(
    IOobject
    (
        "U", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    hf * (Uf & mesh.Sf())
);
U.write();

Info << "Calculating V mass circulation field" << endl;
surfaceScalarField V
(
    IOobject
    (
        "V", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    hf * (Uf & H.delta()),
    VuBCs
);

// fix no flow boundaries
forAll(V.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    const word htype = h.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || utype == "slip" || htype == "geostrophicZeroFlux"
    )
    {
        V.boundaryFieldRef()[patchi] == 0;
    }
    else if(utype == "fixedValue")
    {
        u.boundaryFieldRef()[patchi] == Uf.boundaryField()[patchi];
    }
}
