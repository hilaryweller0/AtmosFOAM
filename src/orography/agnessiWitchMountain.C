#include "agnessiWitchMountain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(agnessiWitchMountain, 0);
addToRunTimeSelectionTable(mountain, agnessiWitchMountain, dict);

agnessiWitchMountain::agnessiWitchMountain(const dictionary& dict)
:
    h0_(dict.lookup("peakHeight")),
    centre_(dict.lookup("mountainCentre")),
    halfwidth_(dict.lookup("mountainHalfWidth"))
{};

dimensionedScalar agnessiWitchMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    
    scalar r = horizontalDist(p, centre_);
    h = h0_/(1 + magSqr(r/halfwidth_.value()));

    return h;
}
}
