#include "WitchOfAgnesi3DMountain.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(WitchOfAgnesi3DMountain, 0);
addToRunTimeSelectionTable(Mountain, WitchOfAgnesi3DMountain, dict);

WitchOfAgnesi3DMountain::WitchOfAgnesi3DMountain(const IOdictionary& dict) :
    a(readScalar(dict.lookup("mountainHalfWidth"))),
    h0(readScalar(dict.lookup("mountainPeakHeight")))
{}


scalar WitchOfAgnesi3DMountain::heightAt(const point& p) const
{
    return h0/Foam::pow(1 + sqr(p.x()/a) + sqr(p.y()/a), 1.5);
}

scalar WitchOfAgnesi3DMountain::gradientAt(const scalar x) const
{
    return 0;
}

scalar WitchOfAgnesi3DMountain::timeToCross(const scalar u0, const scalar H) const
{
    return 0;
}

scalar WitchOfAgnesi3DMountain::halfWidth() const
{
    return a;
}

