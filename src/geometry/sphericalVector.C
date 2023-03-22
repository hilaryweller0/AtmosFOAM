#include "sphericalVector.H"

#include "vector.H"
#include "tensor.H"

Foam::sphericalVector::sphericalVector(const vector v) : v(v) {};

Foam::sphericalVector::sphericalVector
(
        const scalar u,
        const scalar v,
        const scalar w
)
:
    v(vector(u, v, w))
{};

Foam::vector Foam::sphericalVector::toCartesian(const sphericalVector& point) const
{
    if (mag(point.v) < VSMALL) return vector(0,0,0);
    return point.unitTensor().inv() & v;
}

Foam::tensor Foam::sphericalVector::unitTensor() const
{
    const vector kHat(0,0,1);
    const vector rHat = v / mag(v);

    vector latHat;
    if (mag(kHat ^ rHat) < VSMALL)
    {
        if ((kHat & rHat) > VSMALL)
        {
            // North Pole
            latHat = vector(1, 0, 0);
        }
        else
        {
            // South Pole
            latHat = vector(-1, 0, 0);
        }
    }
    else
    {
        latHat = rHat ^ (kHat ^ rHat);
        latHat /= mag(latHat);
    }

    const vector lonHat = latHat ^ rHat;
    // TODO: use tmp<tensor>?
    return tensor(lonHat, latHat, rHat);
}

Foam::sphericalVector Foam::convertToLocal
(
    const point& pc,
    const vector& v
)
{
    sphericalVector p(pc);
    return p.unitTensor() & v;
}

