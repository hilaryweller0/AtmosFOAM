#include "schaerCosMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerCosMountain, 0);
addToRunTimeSelectionTable(mountain, schaerCosMountain, dict);

schaerCosMountain::schaerCosMountain(const dictionary& dict)
:
h0("peakHeight", dimLength, readScalar(dict.lookup("peakHeight"))),
a("a", dimLength, readScalar(dict.lookup("halfWidth"))),
lambda("lambda", dimLength, readScalar(dict.lookup("wavelength"))),
xOffset("xOffset", dimLength, dict.lookupOrDefault("xOffset", 0))
{};

dimensionedScalar schaerCosMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    dimensionedScalar X("X", dimLength, scalar(p.x()));
    if (p.x() >= (xOffset - a).value() && p.x() <= (xOffset + a).value())
    {
        dimensionedScalar x = X - xOffset;
        h = h0 * sqr(Foam::cos(M_PI * x / lambda)) * sqr(Foam::cos(0.5*M_PI*x/a));
    }
    return h;
}

dimensionedScalar schaerCosMountain::start() const
{
    return -a;
}

dimensionedScalar schaerCosMountain::end() const
{
    return a;
}

dimensionedScalar schaerCosMountain::timeToCross
(
    const dimensionedScalar u0,
    const dimensionedScalar H
) const
{
    return timeToCrossIntegralAt(u0, H, end()) - timeToCrossIntegralAt(u0, H, start());
}

dimensionedScalar schaerCosMountain::timeToCrossIntegralAt
(
    const dimensionedScalar u0,
    const dimensionedScalar H,
    const dimensionedScalar x
) const
{
    const dimensionedScalar alpha = M_PI/lambda;
    const dimensionedScalar beta = M_PI/(2*a);

    const dimensionedScalar A = 0.25*
    (
        Foam::sin(2.0 * (alpha + beta) * x) / (alpha + beta) +
        Foam::sin(2.0 * (alpha - beta) * x) / (alpha - beta)
    );
    const dimensionedScalar B = 0.5 *
    (
        Foam::sin(2.0*alpha*x)/alpha +
        Foam::sin(2.0*beta*x)/beta
    );

    return x / u0 - h0/(4.0 * u0 * H) * (x + A + B);
}
