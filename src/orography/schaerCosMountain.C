#include "schaerCosMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerCosMountain, 0);
addToRunTimeSelectionTable(mountain, schaerCosMountain, dict);

schaerCosMountain::schaerCosMountain(const dictionary& dict)
:
h0("peakHeight", dimLength, readScalar(dict.lookup("peakHeight"))),
a(readScalar(dict.lookup("halfWidth"))),
lambda(readScalar(dict.lookup("wavelength")))
{};

dimensionedScalar schaerCosMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    if (mag(p.x()) <= a)
    {
        h = h0 * sqr(Foam::cos(M_PI * p.x() / lambda)) * sqr(Foam::cos(0.5*M_PI*p.x()/a));
    }
    return h;
}
