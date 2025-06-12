#ifndef SINE_WAVE_MOUNTAIN_H
#define SINE_WAVE_MOUNTAIN_H

#include "crossableMountain.H"

namespace Foam
{

class sineWaveMountain final : public crossableMountain
{
    public:
    TypeName("sineWave");
    sineWaveMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;
    dimensionedScalar start() const;
    dimensionedScalar end() const;
    dimensionedScalar timeToCross
    (
        const dimensionedScalar u0, 
        const dimensionedScalar H
    ) const;

    private:
    const dimensionedScalar amplitude_;
    const dimensionedScalar waveLength_;
    const dimensionedScalar xStart_;
    const dimensionedScalar xEnd_;
};

}

#endif

