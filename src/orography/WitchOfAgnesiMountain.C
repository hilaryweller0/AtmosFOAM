#include "WitchOfAgnesiMountain.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(WitchOfAgnesiMountain, 0);
addToRunTimeSelectionTable(Mountain, WitchOfAgnesiMountain, dict);

WitchOfAgnesiMountain::WitchOfAgnesiMountain(const IOdictionary& dict) :
    a(readScalar(dict.lookup("mountainHalfWidth"))),
    h0(readScalar(dict.lookup("mountainPeakHeight")))
{}


scalar WitchOfAgnesiMountain::heightAt(const point& p) const
{
    return h0/Foam::pow(1 + sqr(p.x()/a) + sqr(p.y()/a), 1.5);
}

scalar WitchOfAgnesiMountain::gradientAt(const scalar x) const
{
    return 0;
}

scalar WitchOfAgnesiMountain::timeToCross(const scalar u0, const scalar H) const
{
    return 0;
}

scalar WitchOfAgnesiMountain::halfWidth() const
{
    return a;
}

