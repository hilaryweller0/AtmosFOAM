#ifndef TERRAIN_FOLLOWING_TRANSFORM_H
#define TERRAIN_FOLLOWING_TRANSFORM_H

#include "fvCFD.H"

namespace Foam
{

class terrainFollowingTransform
{
    public:
    declareRunTimeSelectionTable
    (
        autoPtr,
        terrainFollowingTransform,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<terrainFollowingTransform> New(const dictionary& dict);
    virtual point physicalToComputational(const point& p) const = 0;
    virtual point computationalToPhysical(const point& p) const = 0;
    virtual ~terrainFollowingTransform() {};
};

}

#endif
