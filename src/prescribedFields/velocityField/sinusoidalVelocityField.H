#ifndef SINUSOIDAL_VELOCITY_FIELD_H
#define SINUSOIDAL_VELOCITY_FIELD_H

/*
A velocity field which is a cosine in each Cartesian direction. Returns:

    u(x,y,z) = (a_x + b_x*cos(2pi*k_x*(x-phi_x)/L_x))
             * (a_y + b_y*cos(2pi*k_y*(y-phi_y)/L_y))
             * (a_z + b_z*cos(2pi*k_z*(z-phi_z)/L_z))
    v(x,y,z) = (c_x + d_x*cos(2pi*l_x*(x-psi_x)/L_x))
             * (c_y + d_y*cos(2pi*l_y*(y-psi_y)/L_y))
             * (c_z + d_z*cos(2pi*l_z*(z-psi_z)/L_z))
    w(x,y,z) = (e_x + f_x*cos(2pi*m_x*(x-chi_x)/L_x))
             * (e_y + f_y*cos(2pi*m_y*(y-chi_y)/L_y))
             * (e_z + f_z*cos(2pi*m_z*(z-chi_z)/L_z))
*/

#include "divergentVelocityField.H"

namespace Foam
{
class sinusoidalVelocityField final : public divergentVelocityField
{
    public:
    TypeName("sinusoidal");
    sinusoidalVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, scalar time) const;

    // The coeffiecients for different terms
    const scalar a_x, a_y, a_z,
                 b_x, b_y, b_z,
                 c_x, c_y, c_z,
                 d_x, d_y, d_z,
                 e_x, e_y, e_z,
                 f_x, f_y, f_z,
                 k_x, k_y, k_z,
                 l_x, l_y, l_z,
                 m_x, m_y, m_z,
                 L_x, L_y, L_z,
                 phi_x, phi_y, phi_z,
                 psi_x, psi_y, psi_z,
                 chi_x, chi_y, chi_z;
};
}
#endif

