#ifndef GEODESIC_HADLEYLIKE_VELOCITY_FIELD_H
#define GEODESIC_HADLEYLIKE_VELOCITY_FIELD_H

// Velocity field for the Hadley-like 3D advection test case of 
/* Dynamical core model intercomparison project:
 Tracer transport test cases
James Kent, Paul A. Ullrich and Christiane Jablonowski*/

#include "geodesicVelocityField.H"

namespace Foam
{
class geodesicHadleyLikeVelocityField final : public geodesicVelocityField
{
    public:
    TypeName("geodesicHadleyLike");
    geodesicHadleyLikeVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, scalar time) const;
    scalar H_, ztop_, u0_, w0_, tau_, rho0_;
    label K_;
};
}
#endif
