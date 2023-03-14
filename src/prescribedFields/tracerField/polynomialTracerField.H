#ifndef POLYNOMIAL_TRACER_FIELD_H
#define POLYNOMIAL_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class polynomialTracerField final : public tracerField
{
    public:
    TypeName("polynomial");
    polynomialTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    // The polynomial coeffiecients for different terms
    const scalar a_const, 
                 a_x, a_y, a_z,
                 a_xx, a_xy, a_xz, a_yy, a_yz, a_zz,
                 a_xxx, a_xxy, a_xxz, a_xyy, a_xyz, a_xzz, 
                        a_yyy, a_yyz, a_yzz, a_zzz;
    // The maximum and minimum possible values of the polynomial
    const scalar polyMax, polyMin;
};

#endif


