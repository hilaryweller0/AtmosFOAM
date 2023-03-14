#include "schaerExpMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerExpMountain, 0);
addToRunTimeSelectionTable(mountain, schaerExpMountain, dict);
addToRunTimeSelectionTable(dualGradeMountain, schaerExpMountain, dict);

schaerExpMountain::schaerExpMountain(const dictionary& dict)
:
h0("peakHeight", dimLength, readScalar(dict.lookup("peakHeight"))),
a(readScalar(dict.lookup("halfWidth"))),
lambda(readScalar(dict.lookup("wavelength")))
{};

dimensionedScalar schaerExpMountain::coarseHeightAt(const point& p) const
{
    // eqn 14 in doi.org/10.1175/MWR-D-10-05046.1
    return 0.5 * h0 * exp(-sqr(p.x() / a));
}

dimensionedScalar schaerExpMountain::fineHeightAt(const point& p) const
{
    // eqn 46 in doi.org/10.1175/1520-0493(2002)130<2459:ANTFVC>2.0.CO;2
    dimensionedScalar h =  h0 * exp(-sqr(p.x() / a)) * sqr(Foam::cos(M_PI*p.x()/lambda));
    return h - coarseHeightAt(p);
}

dimensionedScalar schaerExpMountain::start() const
{
    return dimensionedScalar("start", dimLength, -a);
}

dimensionedScalar schaerExpMountain::end() const
{
    return dimensionedScalar("end", dimLength, a);
}

