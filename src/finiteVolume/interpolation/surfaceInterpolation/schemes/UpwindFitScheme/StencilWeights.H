#ifndef StencilWeights_H
#define StencilWeights_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "Fit.H"

namespace Foam
{
    class StencilWeights
    {
    public:
        StencilWeights(const fvMesh& mesh, const word prefix);

        void fitted(
            const label faceI, 
            const autoPtr<Fit>& fit
        );

        void write();

    private:
        const fvMesh& mesh;
        label debugFaceI;
        autoPtr<volScalarField> stencilWeights;
        autoPtr<surfaceScalarField> polynomialTerms;
    };
}

#endif
