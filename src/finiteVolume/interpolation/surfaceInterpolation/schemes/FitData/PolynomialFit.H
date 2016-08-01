#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "Basis.H"
#include "fitResult.H"
#include "fitWeights.H"
#include "localStencil.H"

#include <Eigen/Core>

using Eigen::MatrixXd;

namespace Foam
{
    template<class Polynomial>
    class PolynomialFit
    {
    public:
        PolynomialFit
        (
            const direction dimensions,
            const scalar minSingularValueThreshold = 0.2
        );
        
        autoPtr<fitResult> fit
        (
            fitCoefficients& coefficients,
            fitWeights& weights,
            const localStencil& stencil
        );
    
    private:
        uint32_t findStable
        (
                fitCoefficients& coefficients,
                const localStencil& stencil,
                fitWeights& weights
        );

        void findFullRankCandidates
        (
                const List<uint32_t>& targetLengthCandidates,
                const localStencil& stencil,
                const label targetLength,
                List<uint32_t>& fullRankCandidates,
                SortableList<scalar>& fullRankMinSingularValues
        );

        void populateCoefficients
        (
                scalarList& coefficients,
                const localStencil& stencil,
                uint32_t polynomial,
                label termCount,
                const scalarList& weights
        );

        void populateMatrix
        (
                MatrixXd& B,
                const localStencil& stencil,
                uint32_t polynomial
        );

        void populateMatrix
        (
                MatrixXd& B,
                const localStencil& stencil,
                uint32_t polynomial,
                const scalarList& weights
        );

        label numberOfSetBits(uint32_t i);

        const direction dimensions;
        const scalar minSingularValueThreshold;
    };
}

#ifdef NoRepository
#   include "PolynomialFit.C"
#endif

#endif

