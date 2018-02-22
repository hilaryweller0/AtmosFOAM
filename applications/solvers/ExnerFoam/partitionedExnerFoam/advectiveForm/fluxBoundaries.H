// fixed flux boundary conditions types
wordList surfaceBoundaryTypes(mesh.boundaryMesh().size(), "calculated");
wordList volBoundaryTypes(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf[0].boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        surfaceBoundaryTypes[patchi] = "fixedValue";
        volBoundaryTypes[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        surfaceBoundaryTypes[patchi] = "fixedValue";
        volBoundaryTypes[patchi] = "fixedValue";
    }
}

// Fix the boundary conditions
for(label ip = 0; ip < nParts; ip++)
{
    Info << "Old boundary types are "   
         << volFlux.ddt()[ip].boundaryField().types() << endl;
    volFlux.ddt()[ip].boundaryFieldRef()
        == GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary
        (
            mesh.boundary(), volFlux.ddt()[ip], surfaceBoundaryTypes
        );

    Info << "New boundary types are "   
         << volFlux.ddt()[ip].boundaryField().types() << endl;

    forAll(surfaceBoundaryTypes, patchi)
    {
        if (surfaceBoundaryTypes[patchi] == "fixedValue")
        {
            volFlux.ddt()[ip].boundaryFieldRef()[patchi] == 0;
        }
    }

    volFlux[ip].write();
    volFlux.ddt()[ip].write();
}
Info << surfaceBoundaryTypes << endl;
