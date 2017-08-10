#include "doubleConeMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(doubleConeMountain, 0);
addToRunTimeSelectionTable(mountain, doubleConeMountain, dict);

doubleConeMountain::doubleConeMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    radius_(dict.lookup("mountainRadius")),
    h0m_(dict.lookup("peakHeightM")),
    centrem_(dict.lookup("mountainCentreM")),
    radiusm_(dict.lookup("mountainRadiusM"))
{};

dimensionedScalar doubleConeMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    
    scalar r = mag(p - centre_);
    scalar rm = mag(p - centrem_);
    if (r < radius_.value())
    {
        h = h0_ * (radius_.value() - r)/radius_.value();

    }
    if (rm < radius_.value())
    {
        h = h0m_ * (radiusm_.value() - rm)/radiusm_.value();
    }
    return h;
}

