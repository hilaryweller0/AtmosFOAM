#include "PolynomialFit.H"
#include "fitCoefficients.H"
#include "SortableList.H"
#include "SVD.H"

template<class Polynomial>
Foam::PolynomialFit<Polynomial>::PolynomialFit
(
    const direction dimensions,
    const scalar minSingularValueThreshold
)
:
    dimensions(dimensions),
    minSingularValueThreshold(minSingularValueThreshold)
{}

template<class Polynomial>
autoPtr<fitResult> Foam::PolynomialFit<Polynomial>::fit
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
uint32_t Foam::PolynomialFit<Polynomial>::findStable
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
        findFullRankCandidates
        (
                targetLengthCandidates,
                stencil,
                targetLength,
                fullRankCandidates,
                fullRankMinSingularValues
        );
        const labelList& fullRankIndices = fullRankMinSingularValues.indices();

        forAll(fullRankIndices, candidateI)
        {
            uint32_t candidate = fullRankCandidates[fullRankIndices[candidateI]];
            scalarList w(stencil.size(), scalar(1));
            w[0] = 1024;
            w[1] = 1024;

            do
            {
                scalarList coeffs(stencil.size(), scalar(0));
                scalarRectangularMatrix Binv = populateCoefficients(coeffs, stencil, candidate, targetLength, w);
                
                if (stable(coeffs))
                {
//                    applyHighOrderCorrection(coeffs, stencil, w, Binv, candidate);
                    coefficients.copyFrom(coeffs);
                    weights.copyFrom(w);
                    return candidate;
                }

                w[1] /= 2;
            } while (w[1] >= 1);
        }

        targetLength--;
    } while (targetLength > 0);

    return 0;
}

template<class Polynomial>
void PolynomialFit<Polynomial>::findFullRankCandidates
(
        const List<uint32_t>& targetLengthCandidates,
        const localStencil& stencil,
        const label targetLength,
        List<uint32_t>& fullRankCandidates,
        SortableList<scalar>& fullRankMinSingularValues
)
{
    forAll(targetLengthCandidates, candidateI)
    {
        uint32_t candidate = targetLengthCandidates[candidateI];
        scalarRectangularMatrix B(stencil.size(), targetLength);
        populateMatrix(B, stencil, candidate);

        const SVD svd(B);
        const scalar singularValueRatio = svd.S()[findMax(svd.S())] < VSMALL ? 0 : svd.S()[findMin(svd.S())] / svd.S()[findMax(svd.S())];
        if (singularValueRatio >= minSingularValueThreshold)
        {
            fullRankCandidates.append(candidate);
            fullRankMinSingularValues.append(svd.minNonZeroS());
        }
    }

    fullRankMinSingularValues.reverseSort();
}

template<class Polynomial>
bool PolynomialFit<Polynomial>::stable(const scalarList& coefficients)
{
    scalar maxMagP = 0;
    for (int i=2; i < coefficients.size(); i++)
    {
        if (mag(coefficients[i]) > maxMagP) maxMagP = mag(coefficients[i]);
    }

    return  coefficients[0] >= 0.5 && coefficients[0] <= 1 &&
            coefficients[1] > -SMALL && coefficients[1] <= 0.5 &&
            coefficients[0] - coefficients[1] >= maxMagP;
}

template<class Polynomial>
scalarRectangularMatrix PolynomialFit<Polynomial>::populateCoefficients
(
        scalarList& coefficients,
        const localStencil& stencil,
        uint32_t polynomial,
        label termCount,
        const scalarList& weights
)
{
    scalarRectangularMatrix B(stencil.size(), termCount);
    populateMatrix(B, stencil, polynomial, weights);
    const SVD svd(B);
    const scalarRectangularMatrix& Binv = svd.VSinvUt();

    forAll(coefficients, i)
    {
        coefficients[i] = weights[i]*Binv(0,i);
    }

//    applyHighOrderCorrection(coefficients, stencil, weights, Binv, polynomial);

    return Binv;
}

template<class Polynomial>
void PolynomialFit<Polynomial>::applyHighOrderCorrection
(
        scalarList& coefficients,
        const localStencil& stencil,
        const scalarList& weights,
        const scalarRectangularMatrix& Binv,
        const uint32_t terms
)
{
    Info << "\ncorr ";
    forAll(coefficients, i)
    {
        const scalar second_derivative_upwind = weights[i]*Polynomial::secondXderivative(stencil[0], terms, i, Binv);
        const scalar second_derivative_downwind = weights[i]*Polynomial::secondXderivative(stencil[1], terms, i, Binv);
        scalar corr = 1/48.0 * (-3.0 * second_derivative_upwind + second_derivative_downwind);
        Info << corr << " ";
        coefficients[i] += corr;
    }
}

template<class Polynomial>
void PolynomialFit<Polynomial>::populateMatrix
(
        scalarRectangularMatrix& B,
        const localStencil& stencil,
        uint32_t polynomial
)
{
    scalarList weights(stencil.size(), scalar(1));
    populateMatrix(B, stencil, polynomial, weights);
}

template<class Polynomial>
void PolynomialFit<Polynomial>::populateMatrix
(
        scalarRectangularMatrix& B,
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
label Foam::PolynomialFit<Polynomial>::numberOfSetBits(uint32_t i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
