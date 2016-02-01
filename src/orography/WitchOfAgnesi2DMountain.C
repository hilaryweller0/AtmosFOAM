#include "WitchOfAgnesi2DMountain.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(WitchOfAgnesi2DMountain, 0);
addToRunTimeSelectionTable(Mountain, WitchOfAgnesi2DMountain, dict);

WitchOfAgnesi2DMountain::WitchOfAgnesi2DMountain(const IOdictionary& dict) :
    a(readScalar(dict.lookup("mountainHalfWidth"))),
    h0(readScalar(dict.lookup("mountainPeakHeight")))
{}


scalar WitchOfAgnesi2DMountain::heightAt(const point& p) const
{
      return h0/(1.0 + sqr(p.x()/a));
}

scalar WitchOfAgnesi2DMountain::gradientAt(const scalar x) const
{
    return -2.0 * h0 * x/sqr(a*(1+sqr(x/a)));
}

scalar WitchOfAgnesi2DMountain::timeToCross(const scalar u0, const scalar H) const
{
    return 0;
}

scalar WitchOfAgnesi2DMountain::halfWidth() const
{
    return a;
}

