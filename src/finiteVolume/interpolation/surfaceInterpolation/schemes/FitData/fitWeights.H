#ifndef fitWeights_H
#define fitWeights_H

#include "fvCFD.H"

namespace Foam
{
    class fitWeights
    {
    public:
        fitWeights(const fitWeights& source);
        fitWeights(const label cells);

        void setCentralWeight(const scalar weight, const bool pureUpwind);
        void copyFrom(const scalarList& source);

        friend Ostream& operator<< (Ostream& stream, const fitWeights& c);
        scalar& operator[](int i);
        scalar operator[](int i) const;

    private:
        scalarList cellWeights;
    };
}

#endif

