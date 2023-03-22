#ifndef GEODESIC_HADLEYLIKE_VELOCITY_FIELD_H
#define GEODESIC_HADLEYLIKE_VELOCITY_FIELD_H

// Velocity field for the Hadley-like 3D advection test case of 
/* Dynamical core model intercomparison project:
 Tracer transport test cases
James Kent, Paul A. Ullrich and Christiane Jablonowski*/

#include "divergentVelocityField.H"

class geodesicHadleyLikeVelocityField final : public divergentVelocityField
{
    public:
    TypeName("geodesicHadleyLike");
    geodesicHadleyLikeVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, const Time& t) const;
    scalar a_, H_, ztop_, u0_, w0_, tau_;
    label K_;
};

#endif
