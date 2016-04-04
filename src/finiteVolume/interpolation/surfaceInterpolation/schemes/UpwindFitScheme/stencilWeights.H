#ifndef stencilWeights_H
#define stencilWeights_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "fitResult.H"

namespace Foam
{
    class stencilWeights
    {
    public:
        stencilWeights(const fvMesh& mesh, const word prefix);

        void fitted(
            const label faceI, 
            const autoPtr<fitResult>& fit,
            const List<point>& stencil
        );

        void write();

    private:
        void populateStencilWeights
        (
            const fitResult& fit,
            const List<point>& stencil
        );

        void printStencilCoordinates
        (
            const fitResult& fit,
            const List<point>& stencil
        );

        const fvMesh& mesh;
        label debugFaceI;
        autoPtr<volScalarField> weights;
        autoPtr<surfaceScalarField> polynomialTerms;
        autoPtr<surfaceScalarField> badFits;
    };
}

#endif
