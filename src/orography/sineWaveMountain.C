#include "sineWaveMountain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(sineWaveMountain, 0);
addToRunTimeSelectionTable(mountain, sineWaveMountain, dict);

sineWaveMountain::sineWaveMountain(const dictionary& dict)
:
    amplitude_("amplitude", dimLength, readScalar(dict.lookup("amplitude"))),
    waveLength_("waveLength", dimLength, readScalar(dict.lookup("waveLength"))),
    xStart_("xStart", dimLength, readScalar(dict.lookup("xStart"))),
    xEnd_("xEnd", dimLength, readScalar(dict.lookup("xEnd")))
{};

dimensionedScalar sineWaveMountain::heightAt(const point& p) const
{
    dimensionedScalar h("height", dimLength, scalar(0));
    dimensionedScalar X("X", dimLength, scalar(p.x()));
    if (p.x() >= xStart_.value() && p.x() <= xEnd_.value())
    {
        dimensionedScalar x = X - xStart_;
        h = amplitude_ * (1 - Foam::cos(M_PI * x / waveLength_));
    }
    return h;
}

dimensionedScalar sineWaveMountain::start() const
{
    return xStart_;
}

dimensionedScalar sineWaveMountain::end() const
{
    return xEnd_;
}

dimensionedScalar sineWaveMountain::timeToCross
(
    const dimensionedScalar u0,
    const dimensionedScalar H
) const
{
    return dimensionedScalar("", dimTime, scalar(0));
}
}
