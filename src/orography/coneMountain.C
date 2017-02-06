#include "coneMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(coneMountain, 0);
addToRunTimeSelectionTable(mountain, coneMountain, dict);

coneMountain::coneMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    radius_(dict.lookup("mountainRadius"))
{};

dimensionedScalar coneMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    
    scalar r = mag(p - centre_);
    if (r < radius_.value())
    {
        h = h0_ * (radius_.value() - r)/radius_.value();
    }
    return h;
}

