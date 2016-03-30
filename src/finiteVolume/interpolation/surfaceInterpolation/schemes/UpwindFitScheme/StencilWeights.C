#include "StencilWeights.H"
#include "fieldAccess.H"

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

    std::ostringstream weightsFilename;
    weightsFilename << prefix << "Weights" << debugFaceI;

    stencilWeights.reset(new volScalarField
        (
            IOobject
            (
                weightsFilename.str(),
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

    std::ostringstream polynomialTermsFilename;
    polynomialTermsFilename << prefix << "PolynomialTerms";

    polynomialTerms.reset(new surfaceScalarField
        (
            IOobject
            (
                polynomialTermsFilename.str(),
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
        const autoPtr<fitResult>& fit
)
{
    if (faceI == debugFaceI) populateStencilWeights(fit());
    fieldAccess(polynomialTerms(), faceI) = fit->polynomialTerms;
}

void Foam::StencilWeights::write()
{
    if (debugFaceI > -1) stencilWeights->write();
    polynomialTerms->write();
}

void Foam::StencilWeights::populateStencilWeights(const fitResult& fit)
{
    forAll(fit.stencilPoints, stencilI)
    {
        forAll(mesh.C(), cellI)
        {
            if (mesh.C()[cellI] == fit.stencilPoints[stencilI])
            {
                stencilWeights()[cellI] = fit.coeffs[stencilI];
                if (stencilI == 0)
                {
                    stencilWeights()[cellI] += 1.0; // re-add upwind weight
                }
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& boundaryPatch = mesh.boundary()[patchI];
            fvPatchField<scalar>& weightsPatch = stencilWeights->boundaryField()[patchI];
            forAll(boundaryPatch.Cf(), boundaryFaceI)
            {
                if (boundaryPatch.Cf()[boundaryFaceI] == fit.stencilPoints[stencilI])
                {
                    weightsPatch[boundaryFaceI] = fit.coeffs[stencilI];
                    if (stencilI == 0)
                    {
                        weightsPatch[boundaryFaceI] += 1.0; // re-add upwind weight
                    }
                }
            }
        }
    }
}
