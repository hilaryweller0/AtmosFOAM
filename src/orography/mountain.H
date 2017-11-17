#ifndef MOUNTAIN_H
#define MOUNTAIN_H

#include "fvCFD.H"

namespace Foam
{

class mountain
{
    public:
    declareRunTimeSelectionTable
    (
        autoPtr,
        mountain,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<mountain> New(const dictionary& dict);
    virtual dimensionedScalar heightAt(const point& p) const = 0;
    virtual ~mountain() {};
};

inline scalar horizontalDist(const point& p1, const point& p2)
{
    return sqrt(sqr(p1.x() - p2.x()) + sqr(p1.y() - p2.y()));
}

}

#endif

