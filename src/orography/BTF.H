#ifndef BTF_H
#define BTF_H

#include "fvCFD.H"
#include "Mountain.H"

class BTF
{
    public:
    BTF(const Mountain& mountain, const IOdictionary& dict);
    point transform(const point& geometric) const;
    private:
    const Mountain& mountain;
    const scalar H;
};

#endif
