const polyBoundaryMesh& patches = mesh.boundaryMesh();
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
    if (pp.name() == "outlet")
    {
        fvsPatchField<vector>& patch = interpGradThetaf.boundaryFieldRef()[patchI];
        forAll(patch, faceI)
        {
            patch[faceI] = pTraits<vector>::zero;
        }
    }
}

