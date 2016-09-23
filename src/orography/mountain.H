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

}

#endif

