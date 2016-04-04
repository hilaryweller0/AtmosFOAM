#include "fitWeights.H"

Foam::fitWeights::fitWeights
(
    const label cells
)
:
cellWeights(cells, scalar(1)),
constantWeight(1),
xLinearWeight(1)
{}

void Foam::fitWeights::setCentralWeight
(
    const scalar weight,
    const bool pureUpwind
)
{
    cellWeights[0] = weight;
    if (!pureUpwind) cellWeights[1] = weight;
}

void Foam::fitWeights::setConstantWeight(const scalar weight)
{
    constantWeight = weight;
}

void Foam::fitWeights::setXlinearWeight(const scalar weight)
{
    xLinearWeight = weight;
}

void Foam::fitWeights::removeDownwindWeight()
{
    cellWeights[1] = 1.0;
}

scalar Foam::fitWeights::downwind() const
{
    return cellWeights[1];
}

scalar& Foam::fitWeights::downwind()
{
    return cellWeights[1];
}

scalar& Foam::fitWeights::operator[](int i)
{
    return cellWeights[i];
}

scalar Foam::fitWeights::operator[](int i) const
{
    return cellWeights[i];
}

scalar Foam::fitWeights::constant() const
{
    return constantWeight;
}

scalar Foam::fitWeights::xLinear() const
{
    return xLinearWeight;
}
