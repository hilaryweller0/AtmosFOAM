#include "gaussianMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(gaussianMountain, 0);
addToRunTimeSelectionTable(mountain, gaussianMountain, dict);

gaussianMountain::gaussianMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    halfwidth_(dict.lookup("mountainHalfWidth"))
{};

dimensionedScalar gaussianMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    
    scalar r = horizontalDist(p, centre_);
    h = h0_ / pow((1.0 + pow(mag(r)/halfwidth_.value(),2)), 1.5);

    return h;
}
