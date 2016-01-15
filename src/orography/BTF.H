#ifndef BTF_H
#define BTF_H

#include "fvCFD.H"
#include "Mountain.H"

class BTF
{
    public:
    BTF(const autoPtr<Mountain> mountain, const IOdictionary& dict);
    point computationalToPhysical(const point& geometric) const;
    point physicalToComputational(const point& geometric) const;
    private:
    const autoPtr<Mountain> mountain;
    const scalar H;
};

#endif
