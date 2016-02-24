#ifndef TestableUpwindCorrFitData_H
#define TestableUpwindCorrFitData_H

#include "fvCFD.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "FitData.H"

namespace Foam
{

    template<class Polynomial>
    class TestableUpwindCorrFitData
    :
        public FitData
        <
            TestableUpwindCorrFitData<Polynomial>,
            extendedUpwindCellToFaceStencilNew,
            Polynomial
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
            );

            virtual ~TestableUpwindCorrFitData() {}

            void calcFit();

            void calcFit                                                                                      (
                scalarList& coeffsi,
                scalarList& wts,    
                const List<point>&, 
                const scalar wLin,  
                const label faci 
            );
    };

    template<class Polynomial>
    TestableUpwindCorrFitData<Polynomial>::TestableUpwindCorrFitData
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
            TestableUpwindCorrFitData<Polynomial>,
            extendedUpwindCellToFaceStencilNew,
            Polynomial            
        >
        (
            mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
        )
    {}

    template<class Polynomial>
    void TestableUpwindCorrFitData<Polynomial>::calcFit()
    {}

}

#endif
