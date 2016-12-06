#ifndef EARTH_MOUNTAIN_H
#define EARTH_MOUNTAIN_H

#include "fvCFD.H"
#include "mountain.H"
#include "gdal_priv.h"

class earthMountain : public mountain
{
    public:
    TypeName("earth");
    earthMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;
    ~earthMountain();

    private:
    const scalar xResolution;
    const scalar yResolution;
    const direction dimensions;
    const label yIndex;
    GDALDataset* terrain;
    float* scanlines;
    label xSize;
    label ySize;
};

#endif

