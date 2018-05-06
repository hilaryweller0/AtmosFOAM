#include "doubleCosineMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(doubleCosineMountain, 0);
addToRunTimeSelectionTable(mountain, doubleCosineMountain, dict);

doubleCosineMountain::doubleCosineMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    radius_(dict.lookup("mountainRadius")),
    h0m_(dict.lookup("peakHeightM")),
    centrem_(dict.lookup("mountainCentreM")),
    radiusm_(dict.lookup("mountainRadiusM"))
{};

dimensionedScalar doubleCosineMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));

    scalar r = horizontalDist(p, centre_);
    scalar rm = horizontalDist(p, centrem_);
    if (r < radius_.value())
    {
        h = h0_ * 0.5*(1+Foam::cos(M_PI*r/radius_.value()));
    }
    if (rm < radius_.value())
    {
        h = h0m_ * 0.5*(1+Foam::cos(M_PI*rm/radiusm_.value()));
    }
    return h;
}
