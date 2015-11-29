#ifndef BTF_H
#define BTF_H

#include "fvCFD.H"
#include "Mountain.H"

class BTF
{
    public:
    BTF(const Mountain& mountain, const IOdictionary& dict);
    point computationalToPhysical(const point& geometric) const;
    point physicalToComputational(const point& geometric) const;
    private:
    const Mountain& mountain;
    const scalar H;
};

#endif
