#ifndef fitWeights_H
#define fitWeights_H

#include "fvCFD.H"

namespace Foam
{
    class fitWeights
    {
    public:
        fitWeights
        (
            const label cells
        );

        void setCentralWeight(const scalar weight, const bool pureUpwind);
        void setConstantWeight(const scalar weight);
        void setXlinearWeight(const scalar weight);
        void removeDownwindWeight();
        scalar downwind() const;
        scalar& downwind();

        friend Ostream& operator<< (Ostream& stream, const fitWeights& c);
        scalar& operator[](int i);
        scalar operator[](int i) const;
        scalar constant() const;
        scalar xLinear() const;

    private:
        scalarList cellWeights;
        scalar constantWeight;
        scalar xLinearWeight;
    };
}

#endif

