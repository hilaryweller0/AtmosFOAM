#include "schaerExpMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerExpMountain, 0);
addToRunTimeSelectionTable(mountain, schaerExpMountain, dict);

schaerExpMountain::schaerExpMountain(const dictionary& dict)
:
h0("peakHeight", dimLength, readScalar(dict.lookup("peakHeight"))),
a(readScalar(dict.lookup("halfWidth"))),
lambda(readScalar(dict.lookup("wavelength")))
{};

dimensionedScalar schaerExpMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    // eqn 46 in dx.doi.org/10.1175/1520-0493(2002)130<2459:ANTFVC>2.0.CO;2
    return h0 * exp(-sqr(p.x() / a)) * sqr(Foam::cos(M_PI*p.x()/lambda));
}

dimensionedScalar schaerExpMountain::start() const
{
    return dimensionedScalar("start", dimLength, -a);
}

dimensionedScalar schaerExpMountain::end() const
{
    return dimensionedScalar("end", dimLength, a);
}

