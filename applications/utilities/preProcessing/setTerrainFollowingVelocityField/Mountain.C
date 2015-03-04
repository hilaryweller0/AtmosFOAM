#include "Mountain.H"
#include "fvCFD.H"

SchaerCosFineMountain::SchaerCosFineMountain(const scalar lambda) : lambda(lambda) {};

scalar SchaerCosFineMountain::heightAt(const scalar x) const
{
        return sqr(Foam::cos(M_PI * x / lambda));
}

SchaerSmoothMountain::SchaerSmoothMountain(const scalar a, const scalar h0) : a(a), h0(h0) {};

scalar SchaerSmoothMountain::heightAt(const scalar x) const
{
    scalar h = 0;
    if (mag(x) < a)
    {
        h = h0 * sqr(Foam::cos(0.5*M_PI * x / a));
    }
    return h;
}

SchaerCosMountain::SchaerCosMountain(const scalar a, const scalar h0, const scalar lambda) :
    a(a), h0(h0), lambda(lambda)
{
    smooth = new SchaerSmoothMountain(a, h0);
    fine = new SchaerCosFineMountain(lambda);
}

scalar SchaerCosMountain::heightAt(const scalar x) const
{
    return smooth->heightAt(x) * fine->heightAt(x);
}

scalar SchaerCosMountain::gradientAt(const scalar x) const
{
    if (x < -a || x > a)
    {
        return 0.0;
    }
    else
    {
        return - h0 * M_PI * (1/(2*a)*pow(Foam::cos(M_PI * x/lambda), 2) * Foam::sin(M_PI*x/a) +
            pow(Foam::cos(M_PI * x / (2.0*a)), 2) * Foam::sin(2.0*M_PI*x/lambda)/lambda);
    }
}

SchaerCosMountain::~SchaerCosMountain() {
    free(smooth);
    free(fine);
}
