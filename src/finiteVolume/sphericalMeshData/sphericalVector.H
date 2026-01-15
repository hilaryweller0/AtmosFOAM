#ifndef sphericalVector_H
#define sphericalVector_H

#include "point.H"
#include "scalar.H"
#include "tensor.H"
#include "vector.H"

namespace Foam
{

class sphericalVector
{
    public:
    sphericalVector(const vector v);

    sphericalVector
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

    vector toCartesian(const sphericalVector& point) const;
    tensor unitTensor() const;

    //private:
    const vector v;
};

//- convert Cartesian vector to local Cartesian
vector convertToLocal
(
    const point& pc,
    const vector& v
);

}
#endif
