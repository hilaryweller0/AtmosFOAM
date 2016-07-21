#include "stencilWeights.H"
#include "fieldAccess.H"

Foam::stencilWeights::stencilWeights(const fvMesh& mesh, const word prefix)
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

    weights.reset(new volScalarField
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

    std::ostringstream badFitsFilename;
    badFitsFilename << prefix << "BadFits";

    badFits.reset(new surfaceScalarField
        (
            IOobject
            (
                badFitsFilename.str(),
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

    std::ostringstream smallUpwindWeightsFilename;
    smallUpwindWeightsFilename << prefix << "smallUpwindWeights";

    smallUpwindWeights.reset(new surfaceScalarField
        (
            IOobject
            (
                smallUpwindWeightsFilename.str(),
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

    std::ostringstream downwindWeightsFilename;
    downwindWeightsFilename << prefix << "downwindWeights";

    downwindWeights.reset(new surfaceScalarField
        (
            IOobject
            (
                downwindWeightsFilename.str(),
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

void Foam::stencilWeights::fitted(
        const label faceI,
        const autoPtr<fitResult>& fit,
        const List<point>& stencil
)
{
    const point& p = mesh.Cf()[faceI];
//    Info << faceI << " " << p.x() << " " << p.y() << " " << p.z() << " " << 1 + fit->coefficients[0] << " " << fit->coefficients[1] << " " << fit->polynomialTerms << endl;
    
    scalar minP = VGREAT;
    scalar maxP = -VGREAT;
    scalar sumP = 0;
    scalar sumMagP = 0;

    for (int i=2; i < fit->coefficients.size(); i++)
    {
        if (fit->coefficients[i] < minP) minP = fit->coefficients[i];
        if (fit->coefficients[i] > maxP) maxP = fit->coefficients[i];
        sumP += fit->coefficients[i];
        sumMagP += mag(fit->coefficients[i]);
    }
    Info << fit->coefficients[0] << " " << fit->coefficients[1] << " " << minP << " " << maxP << " " << sumP << " " << sumMagP << endl;
    if (faceI == debugFaceI)
    {
        populateStencilWeights(fit(), stencil);
        printStencilCoordinates(fit(), stencil);
        Info << "# coefficients for face " << debugFaceI << " " << fit->coefficients << endl;
        Info << "# polynomialTerms for face " << debugFaceI << " " << fit->polynomialTerms << endl;
        Info << "# cell weights for face " << debugFaceI << " " << fit->weights << endl;
    }
    fieldAccess(polynomialTerms(), faceI) = fit->polynomialTerms;
    fieldAccess(downwindWeights(), faceI) = fit->weights[1];
    fieldAccess(badFits(), faceI) = (fit->good ? 0 : 1);
}

void Foam::stencilWeights::write()
{
    if (debugFaceI > -1) weights->write();
    polynomialTerms->write();
    badFits->write();
    smallUpwindWeights->write();
    downwindWeights->write();
}

void Foam::stencilWeights::populateStencilWeights
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
                weights()[cellI] = fit.coefficients[stencilI];
                if (stencilI == 0)
                {
                    weights()[cellI] += 1.0; // re-add upwind weight
                }
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& boundaryPatch = mesh.boundary()[patchI];
            fvPatchField<scalar>& weightsPatch = weights->boundaryField()[patchI];
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

void Foam::stencilWeights::printStencilCoordinates
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
