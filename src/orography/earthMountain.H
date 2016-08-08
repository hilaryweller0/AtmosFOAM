#ifndef EARTH_MOUNTAIN_MOUNTAIN_H
#define EARTH_MOUNTAIN_MOUNTAIN_H

#include "fvCFD.H"
#include "Mountain.H"
#include "gdal_priv.h"

class earthMountain : public Mountain
{
    public:
    TypeName("earthMountain");
    earthMountain(const IOdictionary& dict);
    scalar heightAt(const point& p) const;
    scalar gradientAt(const scalar x) const;
    scalar timeToCross(const scalar u0, const scalar H) const;
    scalar halfWidth() const;
    ~earthMountain();
    private:
    const scalar xResolution;
    const scalar yResolution;
    const direction dimensions;
    const label yIndex;
    GDALDataset* terrain;
    uint32_t* scanlines;
    label xSize;
    label ySize;
};

#endif

