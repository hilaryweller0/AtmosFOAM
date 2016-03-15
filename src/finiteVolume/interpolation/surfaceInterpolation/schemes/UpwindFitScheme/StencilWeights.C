#include "StencilWeights.H"

Foam::StencilWeights::StencilWeights(const fvMesh& mesh, const word prefix)
:
    mesh(mesh)
{
    debugFaceI = -1;
    mesh.solutionDict()
        .subOrEmptyDict("debug")
        .readIfPresent("interpolationWeightsFaceIndex", debugFaceI, false, false);
    
    // TODO: use this instead of debugFaceI
    //labelList emptyList;
    //labelList debugFaceIlist = debugDict.lookupOrDefault("interpolationWeightsFaceIndices", emptyList, false, false);

    std::ostringstream filename;
    filename << prefix << debugFaceI;

    weightField.reset(new volScalarField
        (
            IOobject
            (
                filename.str(),
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
        mesh,
        0.0,
        "fixedValue"
        )
    );
}

void Foam::StencilWeights::fitted(
        const label faceI,
        const List<scalarList>& coeffs,
        const List<List<point> >& stencilPoints
)
{
    if (faceI != debugFaceI) return;

    forAll(stencilPoints[faceI], stencilI)
    {
        forAll(mesh.C(), cellI)
        {
            if (mesh.C()[cellI] == stencilPoints[faceI][stencilI])
            {
                weightField()[cellI] = coeffs[faceI][stencilI];
                if (stencilI == 0)
                {
                    weightField()[cellI] += 1.0; // re-add upwind weight
                }
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& boundaryPatch = mesh.boundary()[patchI];
            fvPatchField<scalar>& weightsPatch = weightField->boundaryField()[patchI];
            forAll(boundaryPatch.Cf(), boundaryFaceI)
            {
                if (boundaryPatch.Cf()[boundaryFaceI] == stencilPoints[faceI][stencilI])
                {
                    weightsPatch[boundaryFaceI] = coeffs[faceI][stencilI];
                    if (stencilI == 0)
                    {
                        weightsPatch[boundaryFaceI] += 1.0; // re-add upwind weight
                    }
                }
            }
        }
    }
}

void Foam::StencilWeights::write()
{
    weightField->write();
}
