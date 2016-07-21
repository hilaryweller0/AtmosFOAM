#include "PolynomialFit2.H"
#include "MatrixOps.H"
#include "fitCoefficients.H"
#include "SortableList.H"

#include <Eigen/Core>

using Eigen::MatrixXd;

template<class Polynomial>
Foam::PolynomialFit2<Polynomial>::PolynomialFit2
(
    const direction dimensions,
    const scalar minSingularValueThreshold
)
:
    dimensions(dimensions),
    minSingularValueThreshold(minSingularValueThreshold)
{}

template<class Polynomial>
autoPtr<fitResult> Foam::PolynomialFit2<Polynomial>::fit
(
    fitCoefficients& coefficients,
    fitWeights& weights,
    const localStencil& stencil
)
{
    uint32_t polynomial = findStable(coefficients, stencil, weights);
    label terms = numberOfSetBits(polynomial);

    bool goodFit = terms > 0;

    coefficients.applyCorrection(goodFit);

    return autoPtr<fitResult>(new fitResult(
            stencil,
            coefficients,
            weights,
            goodFit,
            polynomial,
            terms
    ));
}

template<class Polynomial>
uint32_t Foam::PolynomialFit2<Polynomial>::findStable
(
        fitCoefficients& coefficients,
        const localStencil& stencil,
        fitWeights& weights
)
{
    label targetLength = min(stencil.size(), Polynomial::nTerms(dimensions));

    List<uint32_t> allCandidates;
    Polynomial::candidates(allCandidates, dimensions);

    do
    {
        List<uint32_t> targetLengthCandidates;
        forAll(allCandidates, candidateI)
        {
            if (numberOfSetBits(allCandidates[candidateI]) == targetLength)
            {
                targetLengthCandidates.append(allCandidates[candidateI]);
            }
        }

        List<uint32_t> fullRankCandidates(0);
        SortableList<scalar> fullRankMinSingularValues(0);

        forAll(targetLengthCandidates, candidateI)
        {
            uint32_t candidate = targetLengthCandidates[candidateI];
	        MatrixXd B(stencil.size(), targetLength);
            populateMatrix(B, stencil, candidate);

            JacobiSVD<MatrixXd> svd(B);
            scalar minSingularValue = svd.singularValues().array().reverse()(0);

            if (minSingularValue >= minSingularValueThreshold)
            {
                fullRankCandidates.append(candidate);
                fullRankMinSingularValues.append(minSingularValue);
            }
        }

        fullRankMinSingularValues.reverseSort();
        const labelList& fullRankIndices = fullRankMinSingularValues.indices();

        forAll(fullRankIndices, candidateI)
        {
            uint32_t candidate = fullRankCandidates[fullRankIndices[candidateI]];
            scalarList w(stencil.size(), scalar(1));
            w[0] = 1000;
            w[1] = 1000;

            do
            {
                scalarList coeffs(stencil.size(), scalar(0));
                populateCoefficients(coeffs, stencil, candidate, targetLength, w);

                scalar maxMagP = 0;
                for (int i=2; i < coeffs.size(); i++)
                {
                    if (mag(coeffs[i]) > maxMagP) maxMagP = mag(coeffs[i]);
                }
                
                if (coeffs[1] < coeffs[0] && coeffs[1] <= 0.5 &&
                    coeffs[0] > 0 && coeffs[1] > -SMALL &&
                    maxMagP*(coeffs.size()-2) < 2)
                {
                    coefficients.copyFrom(coeffs);
                    weights.copyFrom(w);
                    return candidate;
                }

                w[1] -= 1;
            } while (w[1] > 0);
        }

        targetLength--;
    } while (targetLength > 0);

    return 0;
}

template<class Polynomial>
void PolynomialFit2<Polynomial>::populateCoefficients
(
        scalarList& coefficients,
        const localStencil& stencil,
        uint32_t polynomial,
        label termCount,
        const scalarList& weights
)
{
    MatrixXd B(stencil.size(), termCount);
    populateMatrix(B, stencil, polynomial, weights);
    MatrixXd Binv = MatrixOps<MatrixXd>::pseudoInverse(B);

    forAll(coefficients, i)
    {
        coefficients[i] = weights[i]*Binv(0,i);
    }
}

template<class Polynomial>
void PolynomialFit2<Polynomial>::populateMatrix
(
        MatrixXd& B,
        const localStencil& stencil,
        uint32_t polynomial
)
{
    scalarList weights(stencil.size(), scalar(1));
    populateMatrix(B, stencil, polynomial, weights);
}

template<class Polynomial>
void PolynomialFit2<Polynomial>::populateMatrix
(
        MatrixXd& B,
        const localStencil& stencil,
        uint32_t polynomial,
        const scalarList& weights
)
{
    forAll(stencil, stencilI)
    {
        scalar coefficients[Polynomial::nTerms(dimensions)]; 
        Polynomial::addCoeffs(coefficients, stencil[stencilI], 1, dimensions);
        for (label term=0, col=0; term<Polynomial::nTerms(dimensions); term++)
        {
            if ((polynomial & (1 << term)) != 0)
            {
                B(stencilI, col) = coefficients[term] * weights[stencilI];
                col++;
            }
        }
    }
}

template<class Polynomial>
label Foam::PolynomialFit2<Polynomial>::numberOfSetBits(uint32_t i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
