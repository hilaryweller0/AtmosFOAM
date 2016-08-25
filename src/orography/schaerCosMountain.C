#include "schaerCosMountain.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerCosMountain, 0);
addToRunTimeSelectionTable(mountain, schaerCosMountain, dict);

schaerCosMountain::schaerCosMountain(const dictionary& dict)
:
a(readScalar(dict.lookup("halfWidth"))),
h0(readScalar(dict.lookup("peakHeight"))),
lambda(readScalar(dict.lookup("wavelength")))
{};

dimensionedScalar schaerCosMountain::heightAt(const point& p) const
{
    return dimensionedScalar("height", dimLength, scalar(0));
}
