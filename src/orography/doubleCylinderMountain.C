#include "doubleCylinderMountain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(doubleCylinderMountain, 0);
addToRunTimeSelectionTable(mountain, doubleCylinderMountain, dict);

doubleCylinderMountain::doubleCylinderMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    radius_(dict.lookup("mountainRadius")),
    h0m_(dict.lookup("peakHeightM")),
    centrem_(dict.lookup("mountainCentreM")),
    radiusm_(dict.lookup("mountainRadiusM"))
{};

dimensionedScalar doubleCylinderMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));

    scalar r = horizontalDist(p, centre_);
    scalar rm = horizontalDist(p, centrem_);
    if (r < radius_.value())
    {
        h = h0_;
    }
    if (rm < radius_.value())
    {
        h = h0m_;
    }
    return h;
}
}
