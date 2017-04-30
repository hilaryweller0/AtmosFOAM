#include "schaerCosMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerCosMountain, 0);
addToRunTimeSelectionTable(mountain, schaerCosMountain, dict);

schaerCosMountain::schaerCosMountain(const dictionary& dict)
:
h0("peakHeight", dimLength, readScalar(dict.lookup("peakHeight"))),
a(readScalar(dict.lookup("halfWidth"))),
lambda(readScalar(dict.lookup("wavelength"))),
xOffset(dict.lookupOrDefault("xOffset", 0))
{};

dimensionedScalar schaerCosMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    if (p.x() >= xOffset - a && p.x() <= xOffset + a)
    {
        scalar x = p.x() - xOffset;
        h = h0 * sqr(Foam::cos(M_PI * x / lambda)) * sqr(Foam::cos(0.5*M_PI*x/a));
    }
    return h;
}

dimensionedScalar schaerCosMountain::start() const
{
    return dimensionedScalar("start", dimLength, -a);
}

dimensionedScalar schaerCosMountain::end() const
{
    return dimensionedScalar("end", dimLength, a);
}

dimensionedScalar schaerCosMountain::timeToCross
(
    const dimensionedScalar u0,
    const dimensionedScalar H
) const
{
	const dimensionedScalar x("x", dimLength, scalar(2.0 * a));
    const scalar alpha = M_PI/lambda;
    const scalar beta = M_PI/(2.0 * a);

    const dimensionedScalar a("a", dimLength, 0.25*(Foam::sin(2.0 * (alpha + beta) * x.value()) / (alpha + beta) +
         Foam::sin(2.0 * (alpha - beta) * x.value()) / (alpha - beta)));
    const dimensionedScalar b("b", dimLength, 0.5 * (Foam::sin(2.0*alpha*x.value())/alpha + Foam::sin(2.0*beta*x.value())/beta));

    return x / u0 - h0/(4.0 * u0 * H) * (x + a + b);
}
