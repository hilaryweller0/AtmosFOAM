#include "Mountain.H"
#include "fvCFD.H"

SchaerSmoothExpMountain::SchaerSmoothExpMountain(const scalar a, const scalar h0) : a(a), h0(h0) {};

scalar SchaerSmoothExpMountain::heightAt(const scalar x) const
{
    return h0 * Foam::exp(-sqr(x/a));
}

SchaerExpMountain::SchaerExpMountain(const IOdictionary& dict) :
    a(readScalar(dict.lookup("mountainHalfWidth"))),
    h0(readScalar(dict.lookup("mountainPeakHeight"))),
    lambda(readScalar(dict.lookup("mountainWavelength")))
{
    smooth = new SchaerSmoothExpMountain(a, h0);
    fine = new SchaerFineMountain(lambda);
}

SchaerExpMountain::SchaerExpMountain(const scalar a, const scalar h0, const scalar lambda) :
    a(a), h0(h0), lambda(lambda)
{
    smooth = new SchaerSmoothExpMountain(a, h0);
    fine = new SchaerFineMountain(lambda);
}

scalar SchaerExpMountain::heightAt(const scalar x) const
{
    return smooth->heightAt(x) * fine->heightAt(x);
}

scalar SchaerExpMountain::gradientAt(const scalar x) const
{
    return - h0 * Foam::exp(-pow(x/a, 2)) * (M_PI/lambda * Foam::sin(2*M_PI*x/lambda) + 
            2*x / pow(a, 2) * pow(Foam::cos(M_PI*x/lambda), 2));
}

// FIXME: recalculate this
scalar SchaerExpMountain::timeToCross(const scalar u0, const scalar H) const
{
    return 0;
}

scalar SchaerExpMountain::halfWidth() const
{
    return a;
}

SchaerExpMountain::~SchaerExpMountain() {
    free(smooth);
    free(fine);
}
