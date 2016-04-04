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
        const autoPtr<fitResult>& fit,
        const List<point>& stencil
)
{
    if (faceI == debugFaceI)
    {
        populateStencilWeights(fit(), stencil);
        printStencilCoordinates(fit(), stencil);

    }
    fieldAccess(polynomialTerms(), faceI) = fit->polynomialTerms;
}

void Foam::StencilWeights::write()
{
    if (debugFaceI > -1) stencilWeights->write();
    polynomialTerms->write();
}

void Foam::StencilWeights::populateStencilWeights
(
    const fitResult& fit,
    const List<point>& stencil
)
{
    forAll(stencil, stencilI)
    {
        forAll(mesh.C(), cellI)
        {
            if (mesh.C()[cellI] == stencil[stencilI])
            {
                stencilWeights()[cellI] = fit.coefficients[stencilI];
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
                if (boundaryPatch.Cf()[boundaryFaceI] == stencil[stencilI])
                {
                    weightsPatch[boundaryFaceI] = fit.coefficients[stencilI];
                    if (stencilI == 0)
                    {
                        weightsPatch[boundaryFaceI] += 1.0; // re-add upwind weight
                    }
                }
            }
        }
    }
}

void Foam::StencilWeights::printStencilCoordinates
(
    const fitResult& fit,
    const List<point>& stencil
)
{
    Info << "# stencil in mesh coordinates for debugFaceI " << debugFaceI << endl;
    forAll(stencil, i)
    {
        Info << stencil[i].x() << " " << stencil[i].y() << " " << stencil[i].z() << endl;
    }
    Info << endl;

    Info << "# stencil in local coordinates for debugFaceI " << debugFaceI << endl;
    forAll(fit.stencil, i)
    {
        Info << fit.stencil[i].x() << " " << fit.stencil[i].y() << " " << fit.stencil[i].z() << endl;
    }
    Info << endl;
}
