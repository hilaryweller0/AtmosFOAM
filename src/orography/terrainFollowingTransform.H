#ifndef TERRAIN_FOLLOWING_TRANSFORM_H
#define TERRAIN_FOLLOWING_TRANSFORM_H

#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "autoPtr.H"

namespace Foam
{

class terrainFollowingTransform
{
    public:
    TypeName("terrainFollowingTransform");
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
