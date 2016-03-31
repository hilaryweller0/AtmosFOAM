#ifndef stabiliser_H
#define stabiliser_H

#include "weightedMatrix.H"
#include "fitCoefficients.H"

namespace Foam
{
    class stabiliser
    {
    public:
        bool stabilise
        (
            const weightedMatrix& matrix,
            fitCoefficients& coefficients
        ) const;
    };
}

#endif
