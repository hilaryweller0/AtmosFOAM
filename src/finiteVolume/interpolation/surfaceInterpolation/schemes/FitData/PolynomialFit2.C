#include "PolynomialFit2.H"
#include "fitCoefficients.H"
#include "stabiliser.H"

template<class Polynomial>
Foam::PolynomialFit2<Polynomial>::PolynomialFit2
(
    const direction dimensions
)
:
    dimensions(dimensions)
{}

template<class Polynomial>
autoPtr<fitResult> Foam::PolynomialFit2<Polynomial>::fit
(
    fitCoefficients& coefficients,
    fitWeights& weights,
    const localStencil& stencil
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

        Info << "targetLength " << targetLength << endl;
        Info << "targetLengthCandidates " << targetLengthCandidates << endl;

        targetLength--;
    } while (targetLength > 0);

    bool goodFit = (targetLength > 0);

    coefficients.applyCorrection(goodFit);

    return autoPtr<fitResult>(new fitResult(
            stencil,
            coefficients,
            weights,
            goodFit,
            targetLength
    ));
}

template<class Polynomial>
label Foam::PolynomialFit2<Polynomial>::numberOfSetBits(uint32_t i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
