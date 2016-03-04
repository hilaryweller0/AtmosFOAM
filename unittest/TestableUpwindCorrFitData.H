#ifndef TestableUpwindCorrFitData_H
#define TestableUpwindCorrFitData_H

#include "fvCFD.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "cubicUpwindCPCFitPolynomial.H"
#include "FitData.H"

namespace Foam
{
    class TestableUpwindCorrFitData
    :
        public FitData
        <
            TestableUpwindCorrFitData,
            extendedUpwindCellToFaceStencilNew,
            cubicUpwindCPCFitPolynomial
        >
    {
        public:
            TestableUpwindCorrFitData
            (
                const fvMesh& mesh,
                const extendedUpwindCellToFaceStencilNew& stencil,
                const bool linearCorrection,
                const scalar linearLimitFactor,
                const scalar centralWeight
            )
            :
                FitData
                <
                    TestableUpwindCorrFitData,
                    extendedUpwindCellToFaceStencilNew,
                    cubicUpwindCPCFitPolynomial            
                >
                (
                    mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
                )
            {}

            void calcFit() {}

            // TODO: don't know why I have to explicitly delegate to the base class
            // bleugh, polymorphism
            void calcFit
            (
                scalarList& coeffsi,
                scalarList& wts,    
                const List<point>& stencilPoints,
                const scalar wLin,  
                const label faci 
            )
            {
                FitData::calcFit(coeffsi, wts, stencilPoints, wLin, faci);
            }
    };
}

#endif
