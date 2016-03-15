#ifndef StencilWeights_H
#define StencilWeights_H

#include "fvCFD.H"
#include "autoPtr.H"

namespace Foam
{
    class StencilWeights
    {
    public:
        StencilWeights(const fvMesh& mesh, const word prefix);

        void fitted(
            const label faceI, 
            const List<scalarList>& coeffs,
            const List<List<point> >& stencilPoints
        );

        void write();

    private:
        const fvMesh& mesh;
        label debugFaceI;
        autoPtr<volScalarField> weightField;
        //autoPtr<volScalarField> neiWeights;
    };
}

#endif
