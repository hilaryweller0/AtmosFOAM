#ifndef StencilWeights_H
#define StencilWeights_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "fitResult.H"

namespace Foam
{
    class StencilWeights
    {
    public:
        StencilWeights(const fvMesh& mesh, const word prefix);

        void fitted(
            const label faceI, 
            const autoPtr<fitResult>& fit
        );

        void write();

    private:
        void populateStencilWeights(const fitResult& fit);

        const fvMesh& mesh;
        label debugFaceI;
        autoPtr<volScalarField> stencilWeights;
        autoPtr<surfaceScalarField> polynomialTerms;
    };
}

#endif
