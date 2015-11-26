#include "Mountain.H"
#include "fvCFD.H"

SchaerSmoothCosMountain::SchaerSmoothCosMountain(const scalar a, const scalar h0) : a(a), h0(h0) {};

scalar SchaerSmoothCosMountain::heightAt(const point& p) const
{
    scalar h = 0;
    if (mag(p.x()) < a)
    {
        h = h0 * sqr(Foam::cos(0.5*M_PI * p.x() / a));
    }
    return h;
}

SchaerCosMountain::SchaerCosMountain(const IOdictionary& dict) :
    a(readScalar(dict.lookup("mountainHalfWidth"))),
    h0(readScalar(dict.lookup("mountainPeakHeight"))),
    lambda(readScalar(dict.lookup("mountainWavelength")))
{
    smooth = new SchaerSmoothCosMountain(a, h0);
    fine = new SchaerFineMountain(lambda);
}

SchaerCosMountain::SchaerCosMountain(const scalar a, const scalar h0, const scalar lambda) :
    a(a), h0(h0), lambda(lambda)
{
    smooth = new SchaerSmoothCosMountain(a, h0);
    fine = new SchaerFineMountain(lambda);
}

scalar SchaerCosMountain::heightAt(const point& p) const
{
    return smooth->heightAt(p) * fine->heightAt(p);
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

scalar SchaerCosMountain::timeToCross(const scalar u0, const scalar H) const
{
    const scalar x = 2.0 * a;
    const scalar alpha = M_PI/lambda;
    const scalar beta = M_PI/(2.0 * a);

    return x / u0 - h0/(4.0 * u0 * H) * 
        (x + 0.25*(Foam::sin(2.0 * (alpha + beta) * x) / (alpha + beta) +
         Foam::sin(2.0 * (alpha - beta) * x) / (alpha - beta)) +
         0.5 * (Foam::sin(2.0*alpha*x)/alpha + Foam::sin(2.0*beta*x)/beta));
}

scalar SchaerCosMountain::halfWidth() const
{
    return a;
}

SchaerCosMountain::~SchaerCosMountain() {
    free(smooth);
    free(fine);
}
